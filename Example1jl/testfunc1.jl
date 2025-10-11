include("TPSAad.jl")
using .TPSAad

function main()
    # Number of variables
    noi = 5
    
    # Initialize the AD module
    #
    dtpsad_initialize(noi)
    #dtpsad_initialize(noi)
    
    # Initialize variables
    x1 = DTPSAD(); assign_v2!(x1, 3.0, 1)
    x2 = DTPSAD(); assign_v2!(x2, 0.1, 2)
    x3 = DTPSAD(); assign_v2!(x3, 0.1, 3)
    x4 = DTPSAD(); assign_v2!(x4, 1.0, 4)
    x5 = DTPSAD(); assign_v2!(x5, 1.0, 5)
    
    # Test function
    func = (2*cos(x1/x2) + x1/x2 + exp(x2))*x3 + 2*x4 + sinh(x5)
    
    # Alternative test functions (commented out)
    # func = (2*cos(x1/x2) + x1/x2 + exp(x2))*x3
    # func = cos(x1*x2) + x1*x2 + exp(x2)
    
    println("func value and its derivatives w.r.t. 5 variables:")
    println(func.map[1:noi+1])
end

# Run the main function
main()
