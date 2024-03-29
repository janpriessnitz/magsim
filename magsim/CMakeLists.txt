set ( MAGSIM_SOURCES
    Config.cpp
    ConfigReader.cpp
    CoRuCoGenerator.cpp
    Heatbath.cpp
    LatticeGenerator.cpp
    MapReader.cpp
    Metropolis.cpp
    SimulationFactory.cpp
    SpinDynamics.cpp
    SpinLattice.cpp
    Timer.cpp
    TupleReader.cpp
    Util.cpp
)

add_executable(magsim Cli.cpp ${MAGSIM_SOURCES})
target_compile_options(magsim PRIVATE "-O2")
target_compile_options(magsim PRIVATE "-march=native")
target_compile_options(magsim PRIVATE "-fopenmp")

# target_compile_options(magsim PRIVATE "-pg")

target_compile_features(magsim PRIVATE cxx_std_17)


find_package(OpenMP REQUIRED)
target_link_libraries(magsim PRIVATE OpenMP::OpenMP_CXX)
# target_link_options(magsim PRIVATE "-pg")
