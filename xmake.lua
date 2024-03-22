set_project("garm-ls")
set_xmakever("2.7.1")

includes("Common")

add_requires("eigen", "glfw", "glad", "yaml-cpp")

add_rules("mode.debug", "mode.release", "mode.releasedbg")
set_languages("cxxlatest")

option("openmp", {default = true, showmenu = true})
if has_config("openmp") then
    add_requires("openmp")
    add_packages("openmp")
    add_defines("EIGEN_HAS_OPENMP")
end

option("openvdb", {default = false, showmenu = true})
if has_config("openvdb") then
    add_requires("openvdb")
    add_packages("openvdb")
    add_defines("_WITH_OPENVDB")
end

option("disable-rm")
    set_default(false)
    set_showmenu(true)
    add_defines("_DISABLE_RM_")

option("disable-ga")
    set_default(false)
    set_showmenu(true)
    add_defines("_DISABLE_GA_")
    add_deps("disable-rm")
    after_check(function (option)
        if option:dep("disable-rm"):enabled() then
            option:enable(true)
        end
    end)

option("enable-newton")
    set_default(true)
    set_showmenu(true)
    add_defines("_ENABLE_NEWTON_")
    add_deps("disable-ga")
    after_check(function (option)
        if option:dep("disable-ga"):enabled() then
            option:enable(false)
        end
    end)

option("enable-fmm")
    set_default(true)
    set_showmenu(true)
    add_defines("_ENABLE_FMM_")

option("enable-mc")
    set_default(false)
    set_showmenu(true)
    add_defines("_ENABLE_MC_")

option("enable-corr")
    set_default(false)
    set_showmenu(true)
    add_defines("_ENABLE_CORR_")

option("enable-particle")
    set_default(false)
    set_showmenu(true)
    add_defines("_ENABLE_PARTICLE_")

option("use-lerp")
    set_default(false)
    set_showmenu(true)
    add_defines("_USE_LERP_")

option("use-bfecc")
    set_default(false)
    set_showmenu(true)
    add_defines("_USE_BFECC_")

target("core")
    set_kind("static")
    add_files("Source/Core/*.cpp")
    add_includedirs("Source/Core", {public = true})
    add_deps("util")
    add_options("disable-rm", "disable-ga", "enable-newton", "enable-fmm", "enable-mc", "enable-corr", "enable-particle", "use-lerp", "use-bfecc", {public = true})

target("demo")
    set_kind("binary")
    add_files("Source/Demo/*.cpp|standalone.cpp")
    add_deps("core")

target("testadvect2d")
    set_kind("binary")
    set_group("tests")
    add_files("Source/Test/TestAdvect/TestAdvect2D.cpp")
    add_deps("core")

target("testadvect3d")
    set_kind("binary")
    set_group("tests")
    add_files("Source/Test/TestAdvect/TestAdvect3D.cpp")
    add_deps("core")

target("testreinit")
    set_kind("binary")
    set_group("tests")
    add_files("Source/Test/TestReinit/TestReinit.cpp")
    add_deps("core")

target("standalone")
    set_kind("binary")
    add_files("Source/Demo/standalone.cpp")
    add_deps("core")
