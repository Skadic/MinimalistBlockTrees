set(project_EXECUTABLES
        main
        build_bt
        read_bt
        bt_repl)


foreach( exe ${project_EXECUTABLES} )
    add_executable( ${exe} ${exe}.cpp)
endforeach()


foreach( exe ${project_EXECUTABLES} )
    target_link_libraries( ${exe} blocktree )
endforeach()
