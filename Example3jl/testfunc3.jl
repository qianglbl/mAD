include("TPSAad.jl")
using .TPSAad

# Helper function for quadmap
function quadmap!(tpsaPtc::Vector{DTPSAD}, tau::DTPSAD, outin::DTPSAD, xk::Float64)
    tpsaPtcTemp = Vector{DTPSAD}(undef, 4)
    
    # Calculate gamma terms
    gammax = (1.0 + tpsaPtc[2]^2)/tpsaPtc[1]
    gammay = (1.0 + tpsaPtc[4]^2)/tpsaPtc[3]
    
    if xk == 0.0
        tpsaPtc[1] = tpsaPtc[1] - 2*tpsaPtc[2]*tau + tau^2*gammax
        tpsaPtc[3] = tpsaPtc[3] - 2*tpsaPtc[4]*tau + tau^2*gammay
    elseif xk > 0.0
        sqk = sqrt(outin)
        co = cos(sqk*tau)
        si = sin(sqk*tau)
        ch = cosh(sqk*tau)
        sh = sinh(sqk*tau)
        
        tpsaPtcTemp[1] = tpsaPtc[1]*co^2 - 
                         2*tpsaPtc[2]*co*si/sqk + (si/sqk)^2*gammax
        tpsaPtcTemp[2] = tpsaPtc[1]*co*si*sqk + 
                         tpsaPtc[2]*(co^2-si^2) - co*si/sqk*gammax
        tpsaPtcTemp[3] = tpsaPtc[3]*ch^2 - 
                         2*tpsaPtc[4]*ch*sh/sqk + (sh/sqk)^2*gammay
        tpsaPtcTemp[4] = -tpsaPtc[3]*ch*sh*sqk + 
                         tpsaPtc[4]*(ch^2+sh^2) - ch*sh/sqk*gammay
                         
        copyto!(tpsaPtc, tpsaPtcTemp)
    elseif xk < 0.0
        sqk = sqrt(-outin)
        co = cos(sqk*tau)
        si = sin(sqk*tau)
        ch = cosh(sqk*tau)
        sh = sinh(sqk*tau)
        
        tpsaPtcTemp[1] = tpsaPtc[1]*ch^2 - 
                         2*tpsaPtc[2]*ch*sh/sqk + (sh/sqk)^2*gammax
        tpsaPtcTemp[2] = -tpsaPtc[1]*ch*sh*sqk + 
                         tpsaPtc[2]*(ch^2+sh^2) - ch*sh/sqk*gammax
        tpsaPtcTemp[3] = tpsaPtc[3]*co^2 - 
                         2*tpsaPtc[4]*co*si/sqk + (si/sqk)^2*gammay
        tpsaPtcTemp[4] = tpsaPtc[3]*co*si*sqk + 
                         tpsaPtc[4]*(co^2-si^2) - co*si/sqk*gammay
                         
        copyto!(tpsaPtc, tpsaPtcTemp)
    end
end

function main()
    # Number of control variables
    noi = 7
    
    # Initialize the AD module
    dtpsad_initialize(noi)
    
    # Initialize lattice variables
    x1 = DTPSAD(); assign_v2!(x1, 0.2, 1)  # drift length
    x2 = DTPSAD(); assign_v2!(x2, 0.1, 2)  # quad length
    x3 = DTPSAD(); assign_v2!(x3, 29.6, 3) # quad strength
    x4 = DTPSAD(); assign_v2!(x4, 0.4, 4)  # drift length
    x5 = DTPSAD(); assign_v2!(x5, 0.1, 5)  # quad length
    x6 = DTPSAD(); assign_v2!(x6, -29.6, 6) # quad strength
    x7 = DTPSAD(); assign_v2!(x7, 0.2, 7)  # drift length
    
    # Initialize Twiss parameters
    twiss = Vector{DTPSAD}(undef, 4)
    for i in 1:4
        twiss[i] = DTPSAD()
    end
    
    twiss[1].map[1] = 1.0   # beta_x
    twiss[2].map[1] = 0.5   # alpha_x
    twiss[3].map[1] = 1.0   # beta_y
    twiss[4].map[1] = -0.5  # alpha_y
    
    # Track through FODO lattice
    nperd = 1
    for i in 1:nperd
        # drift
        quadmap!(twiss, x1, x2, 0.0)
        # quad
        quadmap!(twiss, x2, x3, x3.map[1])
        # drift
        quadmap!(twiss, x4, x2, 0.0)
        # quad
        quadmap!(twiss, x5, x6, x6.map[1])
        # drift
        quadmap!(twiss, x7, x2, 0.0)
    end
    
    println("Beta_x and its derivatives w.r.t. 7 variables:")
    println(twiss[1].map[1:noi+1])
    println("Alpha_x and its derivatives w.r.t. 7 variables:")
    println(twiss[2].map[1:noi+1])
    println("Beta_y and its derivatives w.r.t. 7 variables:")
    println(twiss[3].map[1:noi+1])
    println("Alpha_y and its derivatives w.r.t. 7 variables:")
    println(twiss[4].map[1:noi+1])
end

# Run the main function
main()
