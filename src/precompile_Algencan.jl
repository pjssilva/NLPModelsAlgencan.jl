function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(Tuple{typeof(Algencan.julia_hlp), Int32, Ptr{Float64}, Int32, Ptr{Float64}, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{UInt8}, Ptr{Int32}})
    precompile(Tuple{typeof(Algencan.julia_gjac), Int32, Ptr{Float64}, Ptr{Float64}, Int32, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Int32, Ptr{UInt8}, Ptr{Int32}})
    precompile(Tuple{typeof(Algencan.find_status), Algencan.AlgencanMathProgModel, Float64, Float64, Float64, Int64})
    precompile(Tuple{typeof(Algencan.julia_hl), Int32, Ptr{Float64}, Int32, Ptr{Float64}, Float64, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Int32, Ptr{UInt8}, Ptr{Int32}})
    precompile(Tuple{typeof(Algencan.treat_lower_bounds), Array{Float64, 1}, Array{Float64, 1}})
    precompile(Tuple{typeof(Algencan.julia_fc), Int32, Ptr{Float64}, Ptr{Float64}, Int32, Ptr{Float64}, Ptr{Int32}})
    precompile(Tuple{typeof(Algencan.option2vparam), Algencan.AlgencanMathProgModel})
end
