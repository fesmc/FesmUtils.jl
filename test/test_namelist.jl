using FesmUtils.Namelist

@testset "Namelist" begin

    dir = mktempdir()
    path(n) = joinpath(dir, n)

    write(path("a.nml"), """
    ! file header
    ! ===== Run setup =====

    &ctrl   ! the control group
        ! group comment 1
        ! group comment 2
        run_step  = "ctrl"     ! first comment
                               ! continued here
        units     = "kg m-2 yr-1" "mm d-1"
        times     = -10, -5, 0, 1
        x         = 1.5d3
        flag      = .true.

        ! about y
        y = 2
        s = "a b"
        e = ""
        ! post comment
    /
    &grp2 ! after the name
        a = 1 2.5
    /

    ! trailing
    ! notes
    """)

    @testset "read" begin
        gs = read_namelist(path("a.nml"))
        @test [g.name for g in gs] == ["ctrl", "grp2"]
        g = gs[1]
        @test g.inline == "the control group"
        @test g.pre == "file header\n===== Run setup ====="
        @test g.comment == "group comment 1\ngroup comment 2"
        @test g.post == "post comment"
        @test [k.name for k in g.keys] == ["run_step", "units", "times", "x", "flag", "y", "s", "e"]
        @test g.keys[1].raw == "\"ctrl\""
        @test g.keys[1].comment == "first comment\ncontinued here"
        @test g.keys[6].pre == "about y"
        @test gs[2].inline == "after the name"
        @test gs[2].after == "trailing\nnotes"
        @test group_dict(gs)["grp2"] === gs[2]
    end

    @testset "write round trip" begin
        gs = read_namelist(path("a.nml"))
        write_namelist(path("b.nml"), gs)
        @test read_namelist(path("b.nml")) == gs
        # once written, the text is a fixed point
        write_namelist(path("c.nml"), read_namelist(path("b.nml")))
        @test read(path("b.nml"), String) == read(path("c.nml"), String)
        # header goes before the first group's comments
        write_namelist(path("d.nml"), gs; header="converted\nby a test")
        hs = read_namelist(path("d.nml"))
        @test hs[1].pre == "converted\nby a test\n" * gs[1].pre
        @test hs[2:end] == gs[2:end]
        # aligned keys and comments
        txt = read(path("b.nml"), String)
        @test occursin("    run_step = \"ctrl\"   ! first comment\n                        ! continued here\n", txt)
        @test_throws ErrorException write_namelist(IOBuffer(), [Group("g", [Key("k", "\"a!b\"")])])
    end

    @testset "values (phase 0 rule)" begin
        @test parse_value("\"kg m-2 yr-1\" \"mm d-1\"") == ["kg m-2 yr-1", "mm d-1"]
        @test parse_value("\"kg m-2 yr-1\"") == "kg m-2 yr-1"
        @test parse_value("'pd' 'lgm'") == ["pd", "lgm"]
        @test parse_value("\"pd\", \"lgm\",") == ["pd", "lgm"]
        @test parse_value("\"\"") == ""
        @test parse_value("") == ""
        @test parse_value("\"\" \"x\"") == ["", "x"]
        @test parse_value("-10, -5, 0, 1") == [-10, -5, 0, 1]
        @test parse_value("0.0080 0.0065") == [0.008, 0.0065]
        @test parse_value("1 2.5") isa Vector{Float64}
        @test parse_value("15e3") === 15000.0
        @test parse_value("1.5d3") === 1500.0
        @test parse_value("-.5") === -0.5
        @test parse_value("10") === 10
        for t in ("T", "True", "true", ".true.", ".TRUE.")
            @test parse_value(t) === true
        end
        for t in ("F", "False", ".false.")
            @test parse_value(t) === false
        end
        @test parse_value("True False") == [true, false]
        @test parse_value("\"a\" 1") == ["a", "1"]          # a quoted element: strings
        @test parse_value("none") == "none"
        @test parse_value("a,b") == "a,b"                   # no blank: one element
        @test parse_value("\"a\",\"b\"") == "a,b"           # as nml: one element, quotes removed
        for v in (1.0, 0.1, 1e-5, -21e3, 3.625e14, 42, true, false, "pd", "kg m-2", ["a", "b c"],
                  [1, 2, 3], [0.008, 0.0065], [true, false])
            @test parse_value(format_value(v)) == v
            @test typeof(parse_value(format_value(v))) == typeof(v)
        end
        @test format_value(["pd", "lgm"]) == "\"pd\" \"lgm\""
        @test format_value(false) == "False"
        @test_throws ErrorException format_value("it's")
    end

    @testset "unsupported lines" begin
        write(path("flat.nml"), "ctrl.dtt = 10.0\n")
        @test_throws ErrorException read_namelist(path("flat.nml"))
        write(path("noeq.nml"), "&g\n    a 1\n/\n")
        @test_throws ErrorException read_namelist(path("noeq.nml"))
        write(path("open.nml"), "&g\n    a = 1\n")
        @test_throws ErrorException read_namelist(path("open.nml"))
        write(path("dot.nml"), "&g\n    a = 1.0   ! x.y = 2\n    f = \"a.b\"\n/\n")
        @test parse_value(read_namelist(path("dot.nml"))[1].keys[2].raw) == "a.b"
    end

    @testset "overrides" begin
        ov = parse_overrides(["sim.relax=False", "forcing.coefs==alpha,dTo", "ctrl.time_end=1000"])
        @test ov[1] == ("sim.relax", ["False"], false)
        @test ov[2] == ("forcing.coefs", ["alpha", "dTo"], true)
        @test_throws ErrorException parse_overrides(["ctrl.dtt=1,2"])
        @test_throws ErrorException parse_overrides(["ctrl.dtt"])
    end

    @testset "resolve" begin
        root = mkpath(joinpath(dir, "repo"))
        mkpath(joinpath(root, "input"))
        mkpath(joinpath(root, ".runme"))
        write(joinpath(root, "input", "defaults.nml"), """
        &forcing
            reference = "pd"    ! the reference snapshot
            coefs     = ""
            terms     = ""
        /
        &coef
            method = "const"
            f_min  = 0.0        ! lower bound
        /
        &term
            coef    = ""
            pattern = "snapshots"
            snap0   = ""
            snap1   = ""
        /
        &snapshot
            tas = ""
            pr  = ""
        /
        &itm
            alb_ice = 0.4
            alb_snow = 0.8
        /
        &host
            model = "itm"
            nml_itm = "itm"
        /
        """)
        write(joinpath(root, "input", "library.nml"), """
        &coef_alpha
            method = "series"
        /
        &term_lgm
            coef  = "alpha"
            snap0 = "pi"
            snap1 = "lgm"
        /
        &term_flat
            pattern = "uniform"
            snap1   = "never"
        /
        &snapshot_lgm
            tas = "pmip_tas"
        /
        &snapshot_pi
            tas = "pi_tas"
        /
        &coef_dTo
            method = "derived"
        /
        """)
        write(joinpath(root, ".runme", "info.json"),
              """{"par_defaults": ["input/defaults.nml"], "par_library": ["input/library.nml"]}""")
        write(joinpath(root, "base.nml"), """
        &ctrl
            dtt = 10.0
            time_end = 15e3
        /
        &forcing
            coefs = "alpha"
            terms = "lgm" "flat"
        /
        &snapshot_pd
            tas = "obs_tas"
        /
        &host
            nml_itm = "itm_mine"
        /
        &itm_mine
            alb_ice = 0.31
        /
        """)

        cfg = load_config(root)
        @test cfg.library == ["input/library.nml"]
        groups, order = resolve(cfg, joinpath(root, "base.nml"),
                                parse_overrides(["ctrl.time_end=1000", "forcing.coefs==alpha,dTo",
                                                 "coef_alpha.f_min=0.5"]))
        @test order[1:5] == ["ctrl", "forcing", "snapshot_pd", "host", "itm_mine"]
        @test Set(order) == Set(["ctrl", "forcing", "snapshot_pd", "host", "itm_mine", "coef_alpha",
                                 "term_lgm", "term_flat", "snapshot_pi", "snapshot_lgm", "coef_dTo"])
        # values and their sources
        @test value_of(groups["ctrl"], "time_end") === 1000.0
        @test groups["ctrl"].entries["time_end"].source == "override"
        @test groups["ctrl"].entries["dtt"].source == "base"
        @test groups["ctrl"].template === nothing
        f = groups["forcing"]
        @test f.template == ("input/defaults.nml", "forcing")
        @test value_of(f, "reference") == "pd"
        @test f.entries["reference"].source == "default"
        @test f.entries["reference"].comment == "the reference snapshot"
        @test value_of(f, "coefs") == ["alpha", "dTo"]
        @test f.entries["coefs"].raw == "\"alpha\" \"dTo\""
        @test f.order == ["reference", "coefs", "terms"]
        a = groups["coef_alpha"]
        @test a.source == "library"
        @test a.file == joinpath(root, "input/library.nml")
        @test a.entries["method"].source == "library"
        @test a.entries["f_min"].source == "override"
        @test value_of(a, "f_min") === 0.5
        @test groups["coef_dTo"].source == "library"
        # the terms' snapshots, only for snapshot terms
        @test value_of(groups["snapshot_pi"], "tas") == "pi_tas"
        @test !haskey(groups, "snapshot_never")
        # the reference snapshot from the base, defaults filled
        @test groups["snapshot_pd"].source == "base"
        @test value_of(groups["snapshot_pd"], "pr") == ""
        # a group without defaults through the schema of the same name
        @test value_of(groups["host"], "model") == "itm"
        # a group named by a key, with its defaults group from the rule
        rules = vcat(FORCING_REFS, [GroupRef(r"^host$", "nml_{model}"; template=("input/defaults.nml", "itm"))])
        g2, _ = resolve(cfg, joinpath(root, "base.nml"); follow=rules)
        im = g2["itm_mine"]
        @test im.template == ("input/defaults.nml", "itm")
        @test value_of(im, "alb_ice") === 0.31
        @test im.entries["alb_ice"].source == "base"
        @test value_of(im, "alb_snow") === 0.8
        @test im.entries["alb_snow"].source == "default"
        # writable groups
        wg = Group(g2["itm_mine"])
        @test [k.name for k in wg.keys] == ["alb_ice", "alb_snow"]
        # errors
        @test_throws ErrorException resolve(cfg, joinpath(root, "base.nml"), parse_overrides(["ctrl.nokey=1"]))
        @test_throws ErrorException resolve(cfg, joinpath(root, "base.nml"), parse_overrides(["nogroup.k=1"]))
        @test_throws ErrorException resolve(cfg, joinpath(root, "base.nml"), parse_overrides(["ctrl.dtt=abc"]))
        write(joinpath(root, "bad.nml"), "&forcing\n    terms = \"missing\"\n/\n")
        @test_throws ErrorException resolve(Config(root; library=["input/library.nml"]), joinpath(root, "bad.nml"))
        # an explicit Config: only the listed schema
        cfg2 = Config(root; schema=[r"^coef_\w+$" => ("input/defaults.nml", "coef")],
                      library=["input/library.nml"])
        g3, _ = resolve(cfg2, joinpath(root, "base.nml"))
        @test g3["forcing"].template === nothing
        @test g3["coef_alpha"].template == ("input/defaults.nml", "coef")
        @test value_of(g3["forcing"], "coefs") == "alpha"
        @test g3["term_lgm"].template === nothing
        @test haskey(g3, "snapshot_lgm")       # no pattern key: a snapshot term
        @test !haskey(g3, "snapshot_never")
    end
end
