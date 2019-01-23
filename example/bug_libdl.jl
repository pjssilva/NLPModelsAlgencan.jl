# Triggers weird bug in Libdl
# Just compile Algencan with MA57 support and copy libalgencan.so to /tmp and run using julia
import Libdl
const algencan_lib_path = "/tmp/libalgencan.so"
for i = 1:10000
    println("*******************", i)
    println("Before load lib  =============================================================")
    @assert !(algencan_lib_path in Libdl.dllist())
    algencandl = Libdl.dlopen(algencan_lib_path)
    @assert algencan_lib_path in Libdl.dllist()
    algencansym = Libdl.dlsym(algencandl, :c_algencan)
    println("After load lib =============================================================")
    Libdl.dlclose(algencandl)
    @assert !(algencan_lib_path in Libdl.dllist())
    println("Unloaded =========================================================")
end    	
