macro swap(a,b)
    esc(Expr(:(=), Expr(:tuple, a, b), Expr(:tuple, b, a)))
end

"""
MRI coil compression via PCA

Given multiple MRI surface coil images (`idata_orig`), use SVD/PCA to find a smaller number (`ncoil`) of virtual coil images.
"""
function ir_mri_coil_compress(idata_orig; ncoil=1)
    x, y, z = size(idata_orig)
    idata = reshape(idata_orig, x*y, z)
    ~, S, V = svd(idata, full=false);

    Vr = V[:, 1:ncoil]     # [z ncoil] compression matrix with rank = ncoil
    odata = idata * Vr     # [*N ncoil] compressed data
    odata = reshape(odata, x, y, ncoil) # [(N) ncoil]

    return odata, Vr, Diagonal(S)
end

"""
Generate radial sampling pattern

Arguments:
 - `n1, n2, n3`: image dimensions
 - `lines`: number of radial lines in sampling pattern
"""
function strucrand(n1, n2, n3, lines);
    samp = zeros(Int, n1, n2, n3)
    # Create an array of N points between -n1/2 and n1/2
    x = range(-n1/2, n1/2-1, length=n1)
    y = range(-n2/2, n2/2-1, length=length(x))
    for frame in 1:n3
        # Loop to traverse the kspace ; 0 to pi in steps of π/lines -- 
        # Succesive frames rotated by a small random angle (π/line)*rand()
        # Also, shift the cartesian locations such that the center
        # is now at (n1/2,n1/2) {Otherwise the center would be (0,0)}
        coords = [(round(Int,x[i]*sin(α)+n1/2+0.5),round(Int,y[i]*cos(α)+n2/2+0.5),frame)
            for i in 1:length(x), α in range(π/line*rand(), step=π/line, length=lines)]
        # Create the sampling pattern
        samp[CartesianIndex.(coords)] .= 1
    end
    return samp  
end

## Nonuniform FFT
# Julia implementation of the following matlab function: https://github.com/JeffFessler/reproduce-l-s-dynamic-mri/blob/master/operators/getEnufft.m

# Basis transform generator for getEnufft
function E_basis(basis::String;
        M::Int64, nt::Int64,
        fov::Union{Tuple{Int64,Int64},Nothing} = nothing,
        N::Union{Tuple{Int64,Int64},Nothing} = nothing,
        ksp::Union{Array{Complex{Float64}},Nothing} = nothing)
    if basis == "dirac"
        Bi = ones(M, nt)
    elseif basis == "sinc" || basis == "dirac*dx"
        Bi = ones(M, nt) * prod(fov ./ N, dim=1)
    elseif basis == "rect" # rect(x/dx) ⟺ |dx| * sinc(dx * u)
        dx = abs.(fov) ./ N # usual default
        Bi = ones(M, nt)
        for id = 1:size(ksp,2)
            if dx[id] ≠ 0
                Bi .*= sinc(dx[id] * ksp[:,id,:])
            end
        end
    else
        error("unknown basis_type \"$basis\"")
    end
    Bi
end

# Note that this function in not fully tested, it is just the raw transcription of the getEnufft Matlab function to Julia. I left it here for sake of completeness, but it the notebooks, I used a simplified version specialized for the given task.
function getEnufft_original(sense_maps::Array{Complex{Float64},3};
    nx::Int64 = size(sense_maps,1), ny::Int64 = size(sense_maps,2), nt::Int64,
    nc::Int64 = size(sense_maps,3), mask::Union{BitArray{3},Nothing} = nothing,
    samp::Array = [], N::Tuple{Int64,Int64} = (nx,ny), fov::Tuple{Int64,Int64} = (22,22),
    basis::String = "dirac", ksp::Array,
    M::Int64 = size(ksp,1),
    # nufft params
    om::Array = [], wi::Array = [], donufft::Bool = false, Jd::Array = [6,6],
    Kd::Tuple{Int64,Int64} = floor.(Int64, N.*1.5), n_shift::Tuple{Int64,Int64} = N .÷ 2)
    
    out_basistransform = E_basis(basis, M=M, nt=nt, fov=fov, N=N, ksp=ksp)
    basistransform = donufft ? out_basistransform :
        reshape(out_basistransform, N..., nt)
    
    if donufft
        # input ksp, om w/ size [M,d,nt], and wi w/ size[M,nt]
        isempty(om) && error("cl: need ksp,om,wi for nufft")
        # construct st for nufft
        size(ksp,3) ≠ nt && error("cl: double check ksp dimension")
        st_forw = Array{Any}(undef, size(ksp,3))
        st_backw = Array{Any}(undef, size(ksp,3))
        for tt = 1:size(ksp,3)
            st_forw[tt], st_backw[tt],_ =
                nufft_init(om[:,:,tt], N, n_shift=collect(n_shift))
        end
        E = LinearMap{Complex{Float64}}(
            x -> begin
                S = zeros(M,nt,nc)
                x = reshape(x, nx, ny, nt, 1)
                for tt=1:nt
                    tmp = x[:,:,tt,:] .* sense_maps
                    S[:,tt,:] = reshape(st_forw[tt](tmp)./sqrt(prod(N)),M,1,nc)
                end,
                S .* basistransform
            end,
            S -> begin
                x = zeros(Complex{Float64}, nx,ny,nt)
                S = reshape(S, M,nt,nc)
                S = S .* conj.(basistransform)
                wi ≠ [] && (S = S .* wi) # cl: from otazo
                for tt = 1:nt # cl: '/sqrt(prod(a.imSize))' from otazo
                    tmp = reshape(st_backw[tt](S[:,tt,:])/sqrt(prod(N)),nx,ny,nc)
                    x[:,:,tt] = sum(tmp .* conj.(sense_maps),dims=3)
                end
                x
            end,
            M*nt*nc, nx*ny*nt
        )
    else
        E = LinearMap{Complex{Float64}}(
            x -> begin
                x = reshape(x, nx, ny, nt, 1)
                S = x .* reshape(sense_maps,nx,ny,1,nc)
                S = fftshift(ifft(ifftshift(S,1),1),1)*√(nx)
                S = fftshift(ifft(ifftshift(S,2),2),2)*√(ny)
                samp ≠ nothing && (S = S .* samp) # cl: samp mask only when cartesian samp
                S
            end,
            S -> begin
                S = S .* conj.(basistransform)
                wi ≠ [] && (S = S .* wi) # cl: from otazo
                s = fftshift(fft(ifftshift(S,1),1),1)/√(nx)
                s = fftshift(fft(ifftshift(s,2),2),2)/√(ny)
                dropdims(sum(s .* reshape(conj.(sense_maps),nx,ny,1,nc),dims=4),dims=4)
            end,
            M*nt*nc, nx*ny*nt
        )
    end
    
    E, out_basistransform
end


function fftshift!(
        output::AbstractArray,
        input::AbstractArray,
        dims::NTuple{N,Int}) where {N}
    
    @assert input !== output "input and output must be two distinct arrays"
    @assert any(dims .> 0) "dims can contain only positive values!"
    @assert any(dims .<= ndims(input)) "dims cannot contain larger value than ndims(input) (=$(ndims(input)))"
    @assert size(output) == size(input) "input and output must have the same size"
    @assert eltype(output) == eltype(input) "input and output must have the same eltype"
    
    shifts = [dim in dims ? size(input, dim) ÷ 2 : 0 for dim in 1:ndims(input)]
    circshift!(output, input, shifts)
    
end

function ifftshift!(
        output::AbstractArray,
        input::AbstractArray,
        dims::NTuple{N,Int}) where {N}
    
    @assert input !== output "input and output must be two distinct arrays"
    @assert any(dims .> 0) "dims can contain only positive values!"
    @assert any(dims .<= ndims(input)) "dims cannot contain larger value than ndims(input) (=$(ndims(input)))"
    @assert size(output) == size(input) "input and output must have the same size"
    @assert eltype(output) == eltype(input) "input and output must have the same eltype"
    
    shifts = [dim in dims ? size(input, dim) ÷ 2 + size(input, dim) % 2 : 0 for dim in 1:ndims(input)]
    circshift!(output, input, shifts)
    
end

fftshift!(output::AbstractArray, input::AbstractArray, dims::Int) =
    fftshift!(output, input, (dims,))

ifftshift!(output::AbstractArray, input::AbstractArray, dims::Int) =
    ifftshift!(output, input, (dims,))