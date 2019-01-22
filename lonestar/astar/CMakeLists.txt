app(astar astar.cpp)

add_test_scale(small astar "${BASEINPUT}/structured/rome99.gr")
add_test_scale(web astar "${BASEINPUT}/random/r4-2e26.gr" -delta 8)
