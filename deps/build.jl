using BinDeps
@BinDeps.setup

libalgencan = library_dependency("libalgencan")
depspath = BinDeps.depsdir(libalgencan)
algencanpath = joinpath(depspath, "src", "algencan-3.1.1")
provides(Sources,
URI("http://www.ime.usp.br/~egbirgin/tango/sources/algencan-3.1.1.tgz"),
libalgencan, unpackedpath="algencan-3.1.1")
srcpath = joinpath(depspath, "src")

# Check if there is an already compiled library
if "ALGENCAN_LIB_DIR" in keys(ENV)
    provides(Binaries, ENV["ALGENCAN_LIB_DIR"], libalgencan, os = :Unix)
else
    if !("MA57_SOURCE" in keys(ENV))
        # HSL is not present, compile Algencan sources only.
        provides(SimpleBuild,
        (@build_steps begin
            @warn "You are installing Algencan.jl without HSL libraries."
            @warn "This might preclude good performance."
            @warn "If you can, try to use HSL."
            @warn "See details in the installation section at https://github.com/pjssilva/Algencan.jl ."

            # Get Algencan sources and unpack
            GetSources(libalgencan)
            `tar xf downloads/algencan-3.1.1.tgz --directory=$srcpath`

            # Build Algencan
            @build_steps begin
                ChangeDirectory(algencanpath)
                `make CFLAGS="-O3 -fPIC" FFLAGS="-O3 -ffree-form -fPIC"`
            end

            # Create the shared library
            @build_steps begin
                ChangeDirectory(BinDeps.depsdir(libalgencan))
                CreateDirectory("usr")
                CreateDirectory("usr/lib")
            end
            @build_steps begin
                ChangeDirectory(algencanpath)
                if Sys.isapple()
                    `gfortran -shared -o ../../usr/lib/libalgencan.dylib -Wl,-all_load lib/libalgencan.a -Wl,-noall_load -lgfortran -lblas -llapack`
                else
                    `gfortran -shared -o ../../usr/lib/libalgencan.so -Wl,--whole-archive lib/libalgencan.a -Wl,--no-whole-archive -lgfortran -lblas -llapack`
                end
            end
        end), libalgencan, os = :Unix
        )
    else
        # HSL is present, compile METIS, MA67 and Algencan
        metispath = joinpath(BinDeps.depsdir(libalgencan), "src", "metis-4.0.3")
        metistarpath = joinpath(BinDeps.depsdir(libalgencan), "downloads", "metis-4.0.3.tar.gz")
        ENV["METISPATH"] = metispath

        ma57_src = ENV["MA57_SOURCE"]
        ma57path = joinpath(BinDeps.depsdir(libalgencan), "src", "hsl_ma57-5.2.0")
        ENV["MA57PATH"] = ma57path

        # Build Algencan with HSL
        provides(SimpleBuild,
            (@build_steps begin
                ChangeDirectory(BinDeps.depsdir(libalgencan))
                # Get Algencan sources and unpack
                GetSources(libalgencan)
                `tar xf downloads/algencan-3.1.1.tgz --directory=$srcpath`
                # Get Metis sources and unpack
                FileDownloader("http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz", metistarpath)
                `tar xf downloads/metis-4.0.3.tar.gz --directory=$srcpath`

                # Unpack HSL sources
                CreateDirectory(ma57path)
                `tar xf $ma57_src --directory=$srcpath`

                # Build Metis
                @build_steps begin
                    ChangeDirectory(metispath)
                    `make COPTIONS="-fPIC -O3"`
                end

                # Build HSL
                @build_steps begin
                    ChangeDirectory(ma57path)
                    `patch -p1 -i../../patches/patch_ma57.txt`
                    `./configure --with-metis=$metispath/libmetis.a --prefix=$ma57path CFLAGS="-fPIC -O3" FCFLAGS="-fPIC -O3" FFLAGS="-fPIC -O3"`
                    `make`
                    `make install`
                end

                # Build Algencan
                @build_steps begin
                    ChangeDirectory(algencanpath)
                    `patch -p1 -i../../patches/patch_algencan.txt`
                    `make`
                end

                # Create the shared library
                @build_steps begin
                    ChangeDirectory(BinDeps.depsdir(libalgencan))
                    CreateDirectory("usr")
                    CreateDirectory("usr/lib")
                end
                @build_steps begin
                    ChangeDirectory(algencanpath)
                    if Sys.isapple()
                        `gfortran -shared -o ../../usr/lib/libalgencan.dylib -Wl,-all_load lib/libalgencan.a -Wl,-noall_load -lgfortran -lblas -llapack`
                    else
                        `gfortran -shared -o ../../usr/lib/libalgencan.so -Wl,--whole-archive lib/libalgencan.a -Wl,--no-whole-archive -lgfortran -lblas -llapack`
                    end
                end
            end), libalgencan, os = :Unix
        )
    end
end

@BinDeps.install Dict(:libalgencan => :libalgencan)
