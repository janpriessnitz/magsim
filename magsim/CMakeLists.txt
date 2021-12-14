set ( MAGSIM_SOURCES
    Metropolis.cpp
    SpinLattice.cpp
    Util.cpp
    Config.cpp
)

add_executable(magsim Cli.cpp ${MAGSIM_SOURCES})
target_compile_options(magsim PRIVATE "-O2")
target_compile_options(magsim PRIVATE "-march=native")