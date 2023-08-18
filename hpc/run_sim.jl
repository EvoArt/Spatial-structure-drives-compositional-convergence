using Pkg
Pkg.activate(".")

using Distributions
using Random
using Combinatorics
using DataFrames
using CSV,StatsBase

function Amat(s,c,alpha,self)
    """
    Genereate adjacency matrix of size s,
    with connectance c
    and interaction strength alpha
    """
    A = zeros(s, s)

    for i in 2:s
        for j in 1:i-1
            if rand(Uniform(0,1)) < c
                ints = rand(1:3)
                if ints == 1
                    A[i,j] = rand(Normal(0,alpha))
                elseif ints == 2
                    A[j,i] = rand(Normal(0,alpha))
                else
                    A[i,j] = rand(Normal(0,alpha))
                    A[j,i] = rand(Normal(0,alpha))
                end
            end
            if i == j
                A[i,j] = -self
            end
        end
    end
    return A
end

function Board(s)
    """
    Create s*s matrix for simulation
    to take place on
    """
    return zeros(Int, s, s)
end

function intc(x)
    """
    get coordinates for interactions
    """
    intCoords = []
    for i in -1:1
        for j in -1:1
            if i * j != i+j
                push!(intCoords, (i,j))
            end
        end
    end
    return intCoords
end


function fitnessModifier(board, A, intCoords,focal)
    """
    Calculate the change in fitness of focal individual
    due to local interactions
    """
    modifier = 0

    partners = [Int(board[p[1],p[2]])  for p in intCoords]

    for partner in partners
        if partner != 0
        modifier += A[focal,partner] #
    end
    end
    return modifier
end

function setBoard(board, density, nSpecies)
    """
    set out initial species distributions
    """
    w = Weights(exp.(randn(nSpecies)))

    for i in 1:size(board)[1]
        for j in 1:size(board)[1]
        if rand(Uniform(0,1))< density
            board[i,j] = sample(1:nSpecies,w)
        end
    end
    end
    return board
end

function checkSpace(board, intCoords)

    if sum([board[i[1],i[2]] == 0 for i in intCoords])!= 0
        return true
    else
        return false
    end
end

function Birth!(board, intCoords,focal)

    coords = [(i[1],i[2]) for i in intCoords if board[i[1],i[2]] == 0]
    coord = coords[rand(1:length(coords))]
    board[coord[1],coord[2]] = focal
end

function cull!(board, maxDens,locs,Amatrix,baseFit,intCoords)
    """
    density regulation
    """
    Size = size(board)[1]
    density = sum([i != 0 for i in board])# current density
   
    if density > maxDens # if/while density is too high
        n = Int(round(density - maxDens))
        ws =[]
        mvs = []
        for (i,j) in locs
            if board[i,j] > 0
                focal = board[i,j] # species id                
                locality = [(i,j) .+ ic for ic in intCoords]
                locality = locality[[l for l in 1:8 if checkbounds(Bool,board,locality[l][1],locality[l][2])]]        
                fit = (baseFit[focal] + fitnessModifier(board, Amatrix, locality,focal))
                push!(mvs,CartesianIndex((i,j)))
                push!(ws,fit)
            end
        end
        deaths = sample(mvs,Weights(1 .+ maximum(ws) .- ws),n, replace = false)
        for i in 1:n
            board[deaths[i]] = 0
        end
    end
end

function disperse!(board, sections)
    """
    Shuffle the peices on the board
    """
    Size = size(board)[1]
    deme = Size Ã· sections # length of individual demes
    demes = [i for i in 1:deme:Size] # list of starting points for demes
    for i in demes
        for j in demes
            # shuffle the individual demes
            #each deme is the value from the list + the deme length
            shuffle!(view(board,i: i+deme-1,j:j+deme-1))
        end
    end
end

function getBeta(board, sections,nSpecies)
    """
    beta partitions for a givven scale
    """
	mat = zeros(sections^2,nSpecies)
    Size = size(board)[1]
    deme = Size Ã· sections # length of individual demes
    demes = [i for i in 1:deme:Size] # list of starting points for demes
	Count = 0
    for i in demes
        for j in demes
	Count += 1
            # shuffle the individual demes
            #each deme is the value from the list + the deme length
            site = board[i: i+deme-1,j:j+deme-1]
		for loc in site
			if loc >0
				mat[Count,loc] = 1
			end
		end

        end
    end
    return [1,2,3] # remnant from when diversity calculations were done on the fly
end




function sim(c, nSpecies, density, boardSize, dispersal, demes, alpha,bLim,maxDens,baseFit,A,self )
    # set up board and adjacency matrix
    board = Board(boardSize)
    board = setBoard(board, density, nSpecies)
    Amatrix = A
    intCoords = intc(1)
    moves = combinations(hcat(1:boardSize,1:boardSize),2) |> collect

    # initialise counters
    births = 0
    richness = length(unique(board))
    gens = 0
    stable_gens = 0
    myMoves = shuffle(unique(moves))
    cullMoves = copy(myMoves)
    #run the simulation
    while births < bLim
        gens+=1
        stable_gens +=1
	shuffle!(myMoves)
	shuffle!(cullMoves)
        #get randomly ordered vector of board indeces
        for (i,j) in myMoves
            mv = 1
                if board[i,j] != 0 #not empty
                    #evaluate fitness
                    focal = board[i,j] # species id

                    locality = [(i,j) .+ ic for ic in intCoords]

                    locality = locality[[l for l in 1:8 if checkbounds(Bool,board,locality[l][1],locality[l][2])]]
                    if rand(Uniform(1,3)) < (baseFit[focal] + fitnessModifier(board, Amatrix, locality,focal))
                        # then give birth if space available
                        if checkSpace(board, locality) == true
                        Birth!(board,locality,focal)
                        births += 1                        
                        end
                    end
                end

        end
        #dispersal
        if rand(Uniform(0,1)) < dispersal
            disperse!(board, demes)
        end
        #density regulation
        cull!(board,maxDens,cullMoves,Amatrix,baseFit,intCoords)
        # check for extinctions
        if length(unique(board)) < richness # if an extinction has occurred
            richness = length(unique(board)) #update richness
            births = 0 # reset counter
            stable_gens = 0
        end
    end
	specmat = [sum(board .==spec) for spec in 1:nSpecies]'
    # return data!
    return [[std(baseFit) density c nSpecies boardSize dispersal demes^2 alpha bLim  maxDens richness gens self],specmat]
end

function threadSim(c, nSpecies, density, boardSize, dispersal, demes, alpha,bLim,maxDens, baseFit,A,self)
# initialise dictionary for storing data
D = Dict()
print(Threads.nthreads())
@sync begin
    for i in 1:10
        Threads.@spawn begin
    D[i] = sim(c, nSpecies, density, boardSize, dispersal, demes, alpha,bLim,maxDens,baseFit .+randn(length(baseFit)) .*0.01 ,A,self)
        end
end
end
dat = Array{Float64}(undef, 0, 13)
vat = Array{Float64}(undef, 0, 100)
for (K, V) in D
dat = vcat(dat, V[1])
vat = vcat(vat, V[2])
end
return dat,vat
end
# run the project
let dat = Array{Float64}(undef, 0, 13)
    vat = Array{Float64}(undef, 0, 100)
    for nSpecies in  [100]
        for boardSize in [500]
            b = [boardSize Ã· d for d in 2:boardSize if boardSize/d == boardSize Ã· d]
            Nom = vcat([:stable_gens, :density ,:c, :nSpecies, :boardSize, :dispersal, :demes, :alpha, :bLim,  :maxDens, :richness, :gens, :self]...)
            for c in [0.9]
                for basefit in [0.1]
                    for alpha in [1.0]
                        for density in [0.8]
                            for bLim in boardSize^2 .* [20]
                                for demes in [1]
                                    for (i,dispersal) in enumerate(repeat([0,0.1],1))
                                        for self in [0]
                                            for carbon in 1:2
                                                baseFit = rand(truncated(Normal(1.5,basefit),0.0,3.0),nSpecies)
                                                A = Amat(nSpecies,c,alpha,self)                                
                                                D =threadSim(c, nSpecies, density, boardSize, dispersal, demes, alpha, bLim ,density*boardSize^2,baseFit,A,self)
                                                range = (i-1)*20+(carbon-1)*10+1:(i-1)*20+(carbon)*10
                                                dat = vcat(dat, D[1])
                                                vat = vcat(vat, D[2])
                                                df = DataFrame(dat,Nom)
                                                CSV.write("res" * "$(ARGS[1]).csv", df)
                                                CSV.write("spec_mat" * "$(ARGS[1]).csv",DataFrame(vat,:auto))
                                            end
                                        end                                   
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
