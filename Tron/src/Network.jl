mutable struct Network
    NumberOfRoads    :: Int                    
    NumberOfVertices :: Int                    
    NumberOfLeaves   :: Int                    
    graph            :: SimpleDiGraph          # oriented graph
    Junctions        :: Vector{Junction}       # list of network's junctions
    Adj_mat          :: Array{Int64,2}         # adjacency matrix of the graph
    list_leaves      :: Vector{Int64}          # list of nodes with no outgoing or no ingoing road
    ρ_leaves_t       :: Vector{Float64}        # table of prescribed densities at graph's leaves
    ρmax_t           :: Vector{Float64}        # table of maximum values of ρ on each road
    vmax_t           :: Vector{Float64}        # table of maximum values of v on each road
    σ_t              :: Vector{Float64}        # table of flux's maximizers on each road
    Ψ                :: Vector{Float64}        # maximal flux on leaves

    function Network(ne, nv, Nl, ρmax_t, vmax_t)
        Ne             = ne
        Nv             = nv
        Junctions      = Vector{Junction}(undef, nv)
        graph          = SimpleDiGraph(nv)
        Adj_mat        = zeros(Int64, nv, nv)
        list_leaves    = zeros(Int64, Nl)
        ρ_leaves_t     = zeros(Float64, Nl)
        Ψ              = zeros(Float64, Nl)
        σ_t            = ρmax_t / 2
      new(Ne, Nv, Nl, graph, Junctions, Adj_mat, list_leaves, ρ_leaves_t, ρmax_t, vmax_t, σ_t, Ψ)
    end

    function Network(Ne, Nv, Nl, graph, Junctions, Adj_mat, list_leaves, ρ_leaves_t, ρmax_t, vmax_t, σ_t, Ψ)
        new(Ne, Nv, Nl, graph, Junctions, Adj_mat, list_leaves, ρ_leaves_t, ρmax_t, vmax_t, σ_t, Ψ)
    end
end # class Network

# ----- methods of class Network -----

function set_list_leaves(self::Network, v::Vector{Int64})
    self.list_leaves .= v
end

function set_ρ_leaves_t(self::Network, v::Vector{Float64})
    self.ρ_leaves_t .= v
end

function init_Network(self::Network)
    @inbounds for (k,e) in enumerate(edges(self.graph))
        ##### Construct Junctions
        vL = src(e)
        vR = dst(e)
        # check if assigned, create if needed
        if ! isassigned(self.Junctions, vL) self.Junctions[vL] = Junction() end
        if ! isassigned(self.Junctions, vR) self.Junctions[vR] = Junction() end
        #push neighbors
        set_out_data(self.Junctions[vL], k)
        set_in_data(self.Junctions[vR], k)
        # Map edge to junctions
        self.Adj_mat[vL,vR] = k
    end
    @inbounds for (l,v) in enumerate(self.list_leaves)
        if length(inneighbors(self.graph, v)) == 0 # entry leaf
            idx_out = self.Adj_mat[v, outneighbors(self.graph, v)[1]]
            self.Ψ[l] = fmaxR(self, self.ρ_leaves_t[l], idx_out) 
        elseif length(outneighbors(self.graph, v)) == 0 # exit leaf 
            #! pas de bc sur le flux aux feuilles sortantes
            idx_in = self.Adj_mat[inneighbors(self.graph, v)[1], v]
            #* bonne valeur pour threeways_graph
            # self.Ψ[l] = fmaxL(self, self.σ_t[idx_in], idx_in)  
            self.Ψ[l] = fmaxL(self, 0, idx_in)  
        end
    end
    for v in vertices(self.graph)
        self.Junctions[v].γin_max = Vector{Real}(undef, max(1, self.Junctions[v].Nin)) 
        self.Junctions[v].γout_max = Vector{Real}(undef, max(1, self.Junctions[v].Nout)) 
    end
end

function flow(r::T, rm::R, vm::S) where{T,R,S <: Real}
    return vm*r*(1-r/rm)
end

function f(self::Network, r::T, i::Int64) where {T<:Real} #! not type-stable
    rm, vm = self.ρmax_t[i], self.vmax_t[i]
    return flow(r, rm, vm)
    # return vm*r*(1-r/rm)
end



function dflow(r::T, rm::R, vm::S) where{T,R,S <: Real}
    return vm*(1-2r/rm)
end

function df(self,r,i) #! not type-stable
    rm, vm = self.ρmax_t[i], self.vmax_t[i]
    return dflow(r, rm, vm)
end



fmaxL(self,r,i) = f(self, max(self.σ_t[i],r), i)  #! regmin ?
fmaxR(self,r,i) = f(self, min(self.σ_t[i],r), i)  #! regmax ?

function fmax(self,r,i,dir) #! not type-stable
    # check wether we are computing γL_max (dir=-1) or γR_max (dir=1)
    if dir == -1
        return fmaxL(self, r, i) #* function barrier
    end
    return fmaxR(self, r, i)
end

function local_lax(ρL::T, ρR::R, fL::S, fR::M, λ::N) where{T,R,S,M,N <: Real}
        return 0.5*(fL + fR) - 0.5*λ*(ρR - ρL)
end

function Fnum(self, ρL::T, ρR::S, i) where{T,S <: Real} #! not type-stable
    # println("@Fnum()") 
    fL = f(self,ρL,i) 
    fR = f(self,ρR,i) 
    dfL = df(self,ρL,i) 
    dfR = df(self,ρR,i) 
    return local_lax(ρL, ρR, fL, fR, max(abs(dfL), abs(dfR)) )
    # λ = max(abs(dfL), abs(dfR)) 
    # return local_lax(ρL, ρR, fL, fR, λ)
end 