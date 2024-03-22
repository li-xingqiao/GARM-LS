add_requires("eigen", "glfw", "glad", "yaml-cpp")

add_rules("mode.debug", "mode.release")
set_languages("cxxlatest")

target("util")
    set_kind("static")
    add_files("Util/*.cpp")
    add_includedirs("Util", {public = true})
    add_packages("eigen", "yaml-cpp", {public = true})

target("viewer")
    set_kind("binary")
    add_files("Viewer/*.cpp")
    add_packages("glfw", "glad")
    add_deps("util")
