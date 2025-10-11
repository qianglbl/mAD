include("TPSAad.jl")
using .TPSAad

# Helper function for quadmap
function quadmap!(tpsaPtc::Vector{DTPSAD}, tau::DTPSAD, outin::DTPSAD, xk::Float64)
    tpsaPtcTemp = Vector{DTPSAD}(undef, 4)
    
    if xk == 0.0
        tpsaPtc[1] = tpsaPtc[1] + tpsaPtc[2]*tau
        tpsaPtc[3] = tpsaPtc[3] + tpsaPtc[4]*tau
    elseif xk > 0.0
        sqk = sqrt(outin)
        co = cos(sqk*tau)
        si = sin(sqk*tau)
        ch = cosh(sqk*tau)
        sh = sinh(sqk*tau)
        
        tpsaPtcTemp[1] = tpsaPtc[1]*co + tpsaPtc[2]*si/sqk
        tpsaPtcTemp[2] = -tpsaPtc[1]*si*sqk + tpsaPtc[2]*co
        tpsaPtcTemp[3] = tpsaPtc[3]*ch + tpsaPtc[4]*sh/sqk
        tpsaPtcTemp[4] = tpsaPtc[3]*sh*sqk + tpsaPtc[4]*ch
        
        copyto!(tpsaPtc, tpsaPtcTemp)
    elseif xk < 0.0
        sqk = sqrt(-outin)
        co = cos(sqk*tau)
        si = sin(sqk*tau)
        ch = cosh(sqk*tau)
        sh = sinh(sqk*tau)
        
        tpsaPtcTemp[1] = tpsaPtc[1]*ch + tpsaPtc[2]*sh/sqk
        tpsaPtcTemp[2] = tpsaPtc[1]*sh*sqk + tpsaPtc[2]*ch
        tpsaPtcTemp[3] = tpsaPtc[3]*co + tpsaPtc[4]*si/sqk
        tpsaPtcTemp[4] = -tpsaPtc[3]*si*sqk + tpsaPtc[4]*co
        
        copyto!(tpsaPtc, tpsaPtcTemp)
    end
end

function main()
    # Number of control variables
    noi = 7
    
    # Initialize the AD module
    dtpsad_initialize(noi)
    
    # Initialize lattice variables
    x1 = DTPSAD(); assign_v2!(x1, 0.2, 1)   # drift length
    x2 = DTPSAD(); assign_v2!(x2, 0.1, 2)   # quad length
    x3 = DTPSAD(); assign_v2!(x3, 29.6, 3)  # quad strength
    x4 = DTPSAD(); assign_v2!(x4, 0.4, 4)   # drift length
    x5 = DTPSAD(); assign_v2!(x5, 0.1, 5)   # quad length
    x6 = DTPSAD(); assign_v2!(x6, -29.6, 6) # quad strength
    x7 = DTPSAD(); assign_v2!(x7, 0.2, 7)   # drift length
    
    # Initialize particle coordinates (x,px,y,py)
    pt = Vector{DTPSAD}(undef, 4)
    for i in 1:4
        pt[i] = DTPSAD()
    end
    
    # Set initial coordinates
    pt[1].map[1] = 1.0e-3   # x
    pt[2].map[1] = 1.0e-3   # px
    pt[3].map[1] = -1.0e-3  # y
    pt[4].map[1] = -1.0e-3  # py
    
    # Track through FODO lattice
    nperd = 1
    for i in 1:nperd
        # drift
        quadmap!(pt, x1, x2, 0.0)
        # quad
        quadmap!(pt, x2, x3, x3.map[1])
        # drift
        quadmap!(pt, x4, x2, 0.0)
        # quad
        quadmap!(pt, x5, x6, x6.map[1])
        # drift
        quadmap!(pt, x7, x2, 0.0)
    end
    
    println("X coordinate and its derivatives w.r.t. 7 variables:")
    println(pt[1].map[1:noi+1])
    println("Px coordinate and its derivatives w.r.t. 7 variables:")
    println(pt[2].map[1:noi+1])
    println("Y coordinate and its derivatives w.r.t. 7 variables:")
    println(pt[3].map[1:noi+1])
    println("Py coordinate and its derivatives w.r.t. 7 variables:")
    println(pt[4].map[1:noi+1])
end

# Run the main function
main()
