set ( MAGSIM_SOURCES
    SpinLattice.cpp
    Util.cpp
    Config.cpp
    MapReader.cpp
    ConfigReader.cpp
    TupleReader.cpp
    SpinDynamics.cpp
    LatticeGenerator.cpp
)

add_executable(magsim2 Cli.cpp ${MAGSIM_SOURCES})
target_compile_options(magsim2 PRIVATE "-O2")
target_compile_options(magsim2 PRIVATE "-march=native")
target_compile_options(magsim2 PRIVATE "-fopenmp")

# target_compile_options(magsim2 PRIVATE "-pg")

target_compile_features(magsim2 PRIVATE cxx_std_17)


find_package(OpenMP REQUIRED)
target_link_libraries(magsim2 PRIVATE OpenMP::OpenMP_CXX)
# target_link_options(magsim2 PRIVATE "-pg")
