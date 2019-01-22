/*
 * This file belongs to the Galois project, a C++ library for exploiting parallelism.
 * The code is being released under the terms of the 3-Clause BSD License (a
 * copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2018, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 */

#include "galois/Galois.h"
#include "galois/AtomicHelpers.h"
#include "galois/Reduction.h"
#include "galois/PriorityQueue.h"
#include "galois/Timer.h"
#include "galois/Timer.h"
#include "galois/graphs/LCGraph.h"
#include "galois/graphs/TypeTraits.h"
#include "llvm/Support/CommandLine.h"

#include "Lonestar/BoilerPlate.h"
#include "astar.h"

#include <iostream>
#include <fstream>

namespace cll = llvm::cl;

static const char* name = "Single Source Shortest Path";
static const char* desc =
    "Computes the shortest path from a source node to all nodes in a directed "
    "graph using a modified chaotic iteration algorithm";
static const char* url = "single_source_shortest_path";

static cll::opt<std::string>
    filename(cll::Positional, cll::desc("<input graph>"), cll::Required);

static cll::opt<unsigned int>
    startNode("startNode",
              cll::desc("Node to start search from (default value 0)"),
              cll::init(0));
static cll::opt<unsigned int>
    reportNode("reportNode",
               cll::desc("Node to report distance to(default value 1)"),
               cll::init(1));
static cll::opt<unsigned int>
    stepShift("delta",
              cll::desc("Shift value for the deltastep (default value 13)"),
              cll::init(13));

enum Algo {
  deltaStep,
  serDelta,
  dijkstra
};

const char* const ALGO_NAMES[] = {"deltaStep",
                                  "serDelta", "dijkstra"
                                  };

static cll::opt<Algo>
    algo("algo", cll::desc("Choose an algorithm:"),
         cll::values(
                     clEnumVal(deltaStep, "deltaStep"),
                     clEnumVal(serDelta, "serDelta"),
                     clEnumVal(dijkstra, "dijkstra"),
                     clEnumValEnd),
         cll::init(deltaStep));

// typedef galois::graphs::LC_InlineEdge_Graph<std::atomic<unsigned int>,
// uint32_t>::with_no_lockable<true>::type::with_numa_alloc<true>::type Graph;
//! [withnumaalloc]
using Graph = galois::graphs::LC_CSR_Graph<SNode, uint64_t>::
    with_no_lockable<true>::type ::with_numa_alloc<true>::type;
//! [withnumaalloc]
typedef Graph::GraphNode GNode;

typedef UpdateRequestCommon<GNode> UpdateRequest;

struct UpdateRequestIndexer {
 unsigned shift;

 template <typename R>
 unsigned int operator()(const R& req) const {
   unsigned int t = req.gScore >> shift;
   return t;
 }
};
/*
struct UpdateRequestIndexer: public std::unary_function<UpdateRequest, unsigned int> {
  unsigned int operator() (const UpdateRequest& val) const {
    unsigned int t = val.gScore >> stepShift;
    return t;
  }
}; */
constexpr static const bool TRACK_WORK          = false;
constexpr static const unsigned CHUNK_SIZE      = 64u;
constexpr static const ptrdiff_t EDGE_TILE_SIZE = 512;

//using SSSP                 = BFS_SSSP<Graph, uint32_t, true, EDGE_TILE_SIZE>;
typedef uint32_t Dist;
//using ReqPushWrap          = SSSP::ReqPushWrap;
//using OutEdgeRangeFn       = SSSP::OutEdgeRangeFn;


struct OutEdgeRangeFn {
 Graph& graph;
 auto operator()(const GNode& n) const {
   return graph.edges(n, galois::MethodFlag::UNPROTECTED);
 }

 auto operator()(const UpdateRequest& req) const {
   return graph.edges(req.n, galois::MethodFlag::UNPROTECTED);
 }
};

struct Vertex;
struct Adj {
    uint32_t n;
    uint64_t d_cm;
};

struct Vertex {
    double lat, lon;  // in RADIANS
    std::vector<Adj> adj;

    // Ephemeral state (used during search)
    Vertex* prev; // nullptr if not visited (not in closed set)
};
uint64_t *offsets;
uint32_t *neighbors;
uint64_t *edgeData;
uint64_t numNodes;
uint64_t numEdges;
Vertex* sGraph;
const double EarthRadius_cm = 637100000.0;

uint64_t dist(Graph& graph, const GNode srcNode, const GNode dstNode) {
    // Use the haversine formula to compute the great-angle radians
    SNode* src = &graph.getData(srcNode);
    SNode* dst = &graph.getData(dstNode);
    double latS = std::sin(src->lat - dst->lat);
    double lonS = std::sin(src->lon - dst->lon);
    double a = latS*latS + lonS*lonS*std::cos(src->lat)*std::cos(dst->lat);
    double c = 2*std::atan2(std::sqrt(a), std::sqrt(1-a));

    uint64_t d_cm = c*EarthRadius_cm;
    //printf("src: (%0.10f, %.10f), latS: %.12f , lonS: %.12f, a: %.12lf, sqrt(a) %.12f,  c: %f d_cm :%d\n", src->lat, src->lon, latS, lonS, a, std::sqrt(a), c, d_cm);
    //printf("cos %.12f %.12f\n", std::cos(src->lat) , std::cos(dst->lat));
    //printf("latS2 %.12f, lonS2 %.12f, C2 %.12f\n", latS*latS, lonS*lonS, lonS,std::cos(src->lat)*std::cos(dst->lat));
    //pls::info("dst: (%0.10f, %0.10f)", dst->lat, dst->lon);
    //usleep(5000);
    return d_cm;

}

bool done __attribute__((aligned(128)));
uint64_t done_gScore __attribute__((aligned(128)));

void deltaStepAlgo(Graph& graph, GNode source) {

  //! [reducible for self-defined stats]
  galois::GAccumulator<size_t> BadWork; // non-final updates
  //! [reducible for self-defined stats]
  galois::GAccumulator<size_t> WLEmptyWork; // tasks that doesn't result in an update

  namespace gwl = galois::worklists;

  using PSchunk = gwl::PerSocketChunkFIFO<CHUNK_SIZE>;
  using OBIM    = gwl::OrderedByIntegerMetric<UpdateRequestIndexer, PSchunk>;

  uint64_t initDist = dist(graph, source, reportNode);
  UpdateRequest init(source, initDist, 0, -1);

  graph.getData(source) = 0;

  galois::InsertBag<UpdateRequest> initBag;
  initBag.push(init);

  galois::for_each(galois::iterate(initBag),
                   [&](const UpdateRequest& req, auto& ctx) {
                     constexpr galois::MethodFlag flag =
                         galois::MethodFlag::UNPROTECTED;

                     SNode& data = graph.getData(req.n, flag);
                     //printf("dequeue %d %ld %ld %ld \n", req.n, req.gScore, req.fScore, req.parent);
                     if (done) {
                        if (req.gScore >= done_gScore) return;
                     }
                     if (req.n == reportNode) {
                         done = true;
                         done_gScore = req.gScore;
                     }
                     if (data.dist <= req.gScore) {
                       if (TRACK_WORK)
                         WLEmptyWork += 1;
                       return;
                     }
                     if (TRACK_WORK) {
                       //! [per-thread contribution of self-defined stats]
                       if (data.dist != DIST_INFINITY) {
                         BadWork += 1;
                       }
                       //! [per-thread contribution of self-defined stats]
                     }


                     unsigned int v;
                     while (req.gScore < (v=data.dist)) {
                         if ( __sync_bool_compare_and_swap(&data.dist, v, req.gScore)) {
                            for(Graph::edge_iterator
                                ii = graph.edge_begin(req.n),
                                ee = graph.edge_end(req.n);
                                ii != ee; ++ii) {
                                    GNode dst = graph.getEdgeDst(ii);
                                    uint64_t d = graph.getEdgeData(ii);
                                    // NOTE: Because heuristic is monotonic, no need to have/check closed set
                                    //assert(a.d_cm >= dist(v, a.n)); // OK
                                    uint64_t nFScore = req.fScore + d;
                                    // Due to limited precision, there may be off-by-ones; if
                                    // so, inherit the parent's timestamp to avoid going back
                                    // in time. This is safe, it's typically a 1-cm-off calculation
                                    uint64_t nGScore = std::max(req.gScore, nFScore + dist(graph, dst, reportNode));
                                    //printf("\tenqueue %d %ld %ld %d \n",dst, nGScore, nFScore, req.n);
                                    ctx.push(UpdateRequest(dst, nGScore, nFScore, req.n));

                              }

                             data.parent = req.parent;
                             break;
                         }
                     }
                   },
                   galois::wl<OBIM>(UpdateRequestIndexer{stepShift}),
                   galois::no_conflicts(), galois::loopname("astar"));

  if (TRACK_WORK) {
    //! [report self-defined stats]
    galois::runtime::reportStat_Single("astar", "BadWork", BadWork.reduce());
    //! [report self-defined stats]
    galois::runtime::reportStat_Single("astar", "WLEmptyWork",
                                       WLEmptyWork.reduce());
  }
}

void serDeltaAlgo(Graph& graph, const GNode& source) {

  SerialBucketWL<UpdateRequest, UpdateRequestIndexer> wl(UpdateRequestIndexer{stepShift});
  ;
  uint64_t initDist = dist(graph, source, reportNode);
  UpdateRequest init(source, initDist, 0, -1);

  graph.getData(source) = 0;

  wl.push(init);

  size_t iter = 0ul;
  while (!wl.empty()) {

    auto& curr = wl.minBucket();

    while (!curr.empty()) {
      ++iter;
      UpdateRequest req = curr.front();
      curr.pop_front();

      SNode& data = graph.getData(req.n);
      //printf("dequeue %d %ld %ld %ld \n", req.n, req.gScore, req.fScore, req.parent);
      if (done) {
         if (req.gScore >= done_gScore) continue;
      }
      if (req.n == reportNode) {
          done = true;
          done_gScore = req.gScore;
      }
      if (data.dist <= req.gScore) {
        // empty work
        continue;
      }

      for(Graph::edge_iterator
          ii = graph.edge_begin(req.n),
          ee = graph.edge_end(req.n);
          ii != ee; ++ii) {
              GNode dst = graph.getEdgeDst(ii);
              uint64_t d = graph.getEdgeData(ii);
              // NOTE: Because heuristic is monotonic, no need to have/check closed set
              //assert(a.d_cm >= dist(v, a.n)); // OK
              uint64_t nFScore = req.fScore + d;
              // Due to limited precision, there may be off-by-ones; if
              // so, inherit the parent's timestamp to avoid going back
              // in time. This is safe, it's typically a 1-cm-off calculation
              uint64_t nGScore = std::max(req.gScore, nFScore + dist(graph, dst, reportNode));
              //printf("\tenqueue %d %ld %ld %d \n",dst, nGScore, nFScore, req.n);
              wl.push(UpdateRequest(dst, nGScore, nFScore, req.n));

        }
       data.parent = req.parent;
       data.dist = req.gScore;
    }

    wl.goToNextBucket();
  }

  galois::runtime::reportStat_Single("SSSP-Serial-Delta", "Iterations", iter);
}

void dijkstraAlgo(Graph& graph, const GNode& source) {

  using WL = galois::MinHeap<UpdateRequest>;
  uint64_t initDist = dist(graph, source, reportNode);
  UpdateRequest init(source, initDist, 0, -1);

  graph.getData(source) = 0;

  WL wl;
  wl.push(init);

  size_t iter = 0;

  while (!wl.empty()) {

    if (done) break;
    //counter += 1;
    UpdateRequest req = wl.pop();
    SNode& data = graph.getData(req.n);
    //printf("dequeue %d %ld %ld %ld \n", req.n, req.gScore, req.fScore, req.parent);
    if (!data.parent) {
        data.parent = req.parent;
        if (req.n == reportNode) {
            done = true;
        } else {

            for(Graph::edge_iterator
                ii = graph.edge_begin(req.n),
                ee = graph.edge_end(req.n);
                ii != ee; ++ii) {
                    // Count equivalent swarm tasks. It is 2E, rather than V+E
                    // since neighbor queue_vertex task is enqueued regardless of
                    // whether its parent is set or not
                    iter +=2;
                    GNode dst = graph.getEdgeDst(ii);
                    uint64_t d = graph.getEdgeData(ii);
                    // NOTE: Because heuristic is monotonic, no need to have/check closed set
                    SNode& nData = graph.getData(dst);
                    if (!nData.parent) {
                        //assert(a.d_cm >= dist(v, a.n)); // OK
                        uint64_t nFScore = req.fScore + d;
                        // Due to limited precision, there may be off-by-ones; if
                        // so, inherit the parent's timestamp to avoid going back
                        // in time. This is safe, it's typically a 1-cm-off calculation
                        uint64_t nGScore = std::max(req.gScore, nFScore + dist(graph, dst, reportNode));
                        //printf("\tenqueue %d %ld %ld %d \n",dst, nGScore, nFScore, req.n);
                        wl.push(UpdateRequest(dst, nGScore, nFScore, req.n));
                    }

             }
        }
    }
  }

  galois::runtime::reportStat_Single("SSSP-Dijkstra", "Iterations", iter);
}

uint64_t neighDist(Graph& graph, GNode& v, GNode& w) {

    for (Graph::edge_iterator
            ii = graph.edge_begin(v),
            ee = graph.edge_end(v); ii != ee; ++ii) {
        GNode neighbor = graph.getEdgeDst(ii);
        uint64_t d = graph.getEdgeData(ii);
        if (neighbor == w) return d;
    }
    assert(false);  // w should be in v's adjacency list
    return 0; // make compiler happy
}

void LoadGraph(const char* file) {
    const uint32_t MAGIC_NUMBER = 0x150842A7 + 0;  // increment every time you change the file format
    std::ifstream f;
    f.open(file, std::ios::binary);
    if (!f.is_open()) {
        printf("ERROR: Could not open input file\n");
        exit(1);
    }

    auto readU = [&]() -> uint32_t {
        union U {
            uint32_t val;
            char bytes[sizeof(uint32_t)];
        };
        U u;
        f.read(u.bytes, sizeof(uint32_t));
        assert(!f.fail());
        return u.val;
    };


    auto readD = [&]() -> double {
        union U {
            double val;
            char bytes[sizeof(double)];
        };
        U u;
        f.read(u.bytes, sizeof(double));
        assert(!f.fail());
        return u.val;
    };

    uint32_t magic = readU();
    if (magic != MAGIC_NUMBER) {
        printf("ERROR: Wrong input file format (magic number %d, expected %d)\n",
                magic, MAGIC_NUMBER);
        exit(1);
    }

    numNodes = readU();
    printf("Reading %ld nodes...\n", numNodes);

    sGraph = new Vertex[numNodes];
    uint32_t i = 0;
    while (i < numNodes) {
        sGraph[i].lat = readD();
        sGraph[i].lon = readD();
        uint32_t n = readU();
        sGraph[i].adj.resize(n);
        for (uint32_t j = 0; j < n; j++) sGraph[i].adj[j].n = readU();
        for (uint32_t j = 0; j < n; j++) sGraph[i].adj[j].d_cm = readD()*EarthRadius_cm;
        i++;
    }

    f.get();
    assert(f.eof());

#if 1
    // Print graph

    for (uint32_t i = 0; i < 2/*numNodes*/; i++) {
        printf("%6d: %7f %7f", i, sGraph[i].lat, sGraph[i].lon);
        for (auto a: sGraph[i].adj) printf(" %5d %7ld", a.n, a.d_cm);
        printf("\n");
    }
#endif

    uint64_t adjs = 0;
    for (uint32_t i = 0; i < numNodes; i++) adjs += sGraph[i].adj.size();
    printf("Read %ld nodes, %ld adjacencies\n", numNodes, adjs);
    numEdges = adjs;

    offsets = (uint64_t*) malloc(numNodes * sizeof(uint64_t));
    neighbors = (uint32_t*) malloc(numEdges * sizeof(uint32_t));
    edgeData = (uint64_t*) malloc(numEdges * sizeof(uint64_t));

    uint32_t edge_id = 0;
    for (i=0;i<numNodes;i++) {
        for (uint32_t j=0;j<sGraph[i].adj.size();j++) {
            Adj adj = sGraph[i].adj[j];
            neighbors[edge_id] = adj.n;
            edgeData[edge_id] = adj.d_cm;
            edge_id++;
        }
        offsets[i] = edge_id;
    }
}

void initGraph(Graph& graph, GNode& source, GNode& report) {
  //graph.structureFromFile(filename);
  LoadGraph(filename.c_str());
  galois::graphs::FileGraph fg;
  void* gedgeData =
  fg.fromArrays(offsets , numNodes, (void*) neighbors, numEdges, (char*) edgeData, 8, 0, 0, false, 1);
  memcpy(gedgeData, edgeData, 8 * numEdges);
  graph.allocateFrom(fg);
  graph.constructFrom(fg,0,1);
  std::cout << "Source: " << startNode << " Dest: " << reportNode << "\n";

  /*
  source = *(fg.begin()+1);
  galois::graphs::FileGraph::edge_iterator eb = fg.edge_begin(source);
  while(eb != fg.edge_end(source)) {
     printf("%d %d %ld\n", *eb, fg.getEdgeDst(eb), fg.getEdgeData<uint64_t>(eb));
     eb++;
     }
   for (int i=0;i<10;i++)
      printf("edgeData %d %d \n",i, graph.getEdgeData(i));
   */
  done = false;
  done_gScore = ~0;
  unsigned int id = 0;
  bool foundReport = false;
  bool foundSource = false;
  source = *graph.begin();
  report = *graph.begin();
     printf("%d\n", *(graph.end()) );
  for (Graph::iterator src = graph.begin(), ee =
      graph.end(); src != ee; ++src) {
    SNode& node = graph.getData(*src);
    node.dist = DIST_INFINITY;
    node.lat = sGraph[id].lat;
    node.lon = sGraph[id].lon;
    node.id = id++;
    if (node.id == startNode) {
      source = *src;
      foundSource = true;
      if (graph.edge_begin(source) == graph.edge_end(source)) {
	std::cerr << "ERROR: source has no neighbors\n";
        assert(0);
	abort();
      }
    }
    if (node.id == reportNode) {
      foundReport = true;
      report = *src;
    }
  }

  if (!foundReport) {
    std::cerr << "ERROR: failed to set report (" << reportNode << ").\n";
    assert(0);
    abort();
  }

  if (!foundSource) {
    std::cerr << "ERROR: failed to set source (" << startNode << ").\n";
    assert(0);
    abort();
  }
}

int main(int argc, char** argv) {
  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url);

  Graph graph;
  GNode source, report;
  initGraph(graph, source, report);
  std::cout << "Read " << graph.size() << " nodes\n";


  if (startNode >= graph.size() || reportNode >= graph.size()) {
    std::cerr << "failed to set report: " << reportNode
              << " or failed to set source: " << startNode << "\n";
    assert(0);
    abort();
  }

  auto it = graph.begin();
  std::advance(it, startNode);
  source = *it;
  it     = graph.begin();
  std::advance(it, reportNode);
  report = *it;

  size_t approxNodeData = graph.size() * 64;
  // size_t approxEdgeData = graph.sizeEdges() * sizeof(typename
  // Graph::edge_data_type) * 2;
  galois::preAlloc(numThreads +
                   approxNodeData / galois::runtime::pagePoolSize());
  galois::reportPageAlloc("MeminfoPre");

  if (algo == deltaStep || algo == serDelta) {
    std::cout << "INFO: Using delta-step of " << (1 << stepShift) << "\n";
    std::cout
        << "WARNING: Performance varies considerably due to delta parameter.\n";
    std::cout
        << "WARNING: Do not expect the default to be good for your graph.\n";
  }


  graph.getData(source).dist = 0;

  std::cout << "Running " << ALGO_NAMES[algo] << " algorithm" << std::endl;

  galois::StatTimer Tmain;
  Tmain.start();

  switch (algo) {
  case deltaStep:
    deltaStepAlgo(graph, source);
    break;
  case serDelta:
    serDeltaAlgo(graph, source);
    break;
  case dijkstra:
    dijkstraAlgo(graph, source);
    break;
  default:
    std::abort();
  }

  Tmain.stop();

  galois::reportPageAlloc("MeminfoPost");

  std::vector<GNode> path;

  GNode cur = report;
  while (true) {
      path.push_back(cur);
      if (cur == source) break;
      cur = graph.getData(cur).parent;
      assert(cur);
  }
  std::reverse(path.begin(), path.end());
  uint64_t totalDist_cm = 0;
  for (uint32_t i = 0; i < path.size()-1; i++) {
      uint64_t curDist_cm = neighDist(graph, path[i], path[i+1]);
      totalDist_cm += curDist_cm;
      printf("%4d: %9d -> %9d | %8ld.%02ld m | %8ld.%02ld m\n", i,
              graph.getData(path[i]).id,
              graph.getData(path[i+1]).id,
              curDist_cm / 100, curDist_cm % 100, totalDist_cm / 100, totalDist_cm % 100);
  }

  return 0;
}
