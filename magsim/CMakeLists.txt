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
    XyzReader.cpp
)

find_package(OpenMP REQUIRED)

add_executable(magsim_CoRuCo CoRuCo.cpp ${MAGSIM_SOURCES})
target_compile_options(magsim_CoRuCo PRIVATE "-O2")
target_compile_options(magsim_CoRuCo PRIVATE "-march=native")
target_compile_options(magsim_CoRuCo PRIVATE "-fopenmp")
target_compile_features(magsim_CoRuCo PRIVATE cxx_std_17)

target_link_libraries(magsim_CoRuCo PRIVATE OpenMP::OpenMP_CXX)
# target_link_options(magsim_CoRuCo PRIVATE "-pg")


add_executable(magsim_manual Manual.cpp ${MAGSIM_SOURCES})
target_compile_options(magsim_manual PRIVATE "-O2")
target_compile_options(magsim_manual PRIVATE "-march=native")
target_compile_options(magsim_manual PRIVATE "-fopenmp")
target_compile_features(magsim_manual PRIVATE cxx_std_17)

target_link_libraries(magsim_manual PRIVATE OpenMP::OpenMP_CXX)
# target_link_options(magsim_manual PRIVATE "-pg")