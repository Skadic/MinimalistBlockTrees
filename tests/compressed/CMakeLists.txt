set(test_NAMES
        CBlockTreeFixture
        CBitBlockTreeFixture)


foreach( exe ${test_NAMES} )
    add_executable( check_${exe} ${exe}.cpp)
endforeach()


foreach( exe ${test_NAMES} )
    target_link_libraries(check_${exe} GTest::GTest GTest::Main)
    target_link_libraries(check_${exe} blocktree)
endforeach()

foreach( exe ${test_NAMES} )
    gtest_discover_tests(check_${exe})
endforeach()