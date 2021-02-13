using DSP
using FFTW

mutable struct Signal <: AbstractArray{Union{Real,Complex},1}
      n1::Int
      n2::Int
      vals::Vector{T} where T <: Union{Real, Complex}
      sample_rate::Int

      function Signal(n1::Int, n2::Int, vals, sample_rate=1)
            return new(n1,n2,vals,sample_rate)
      end
end

Base.IndexStyle(::Type{<:Signal}) = IndexLinear()

Base.firstindex(A::Signal) = A.vals[1]
Base.lastindex(A::Signal) = A.vals[end]

Base.size(A::Signal) = (A.n2 - A.n1 + 1,)
Base.length(A::Signal) = length(A.vals)

function Base.getindex(A::Signal, n::Int)
      if n in A.n1:A.n2
            return A.vals[n - A.n1 + 1]
      else
            return 0
      end
end   

function Base.setindex!(A::Signal, v, i::Int)
   
      if i <= A.n1
            A.vals = push!(zeros(typeof(A.vals[1]),A.n1 - i),A.vals...)
            A.vals[1] = v
            A.n1 = i
      elseif A.n1 <= i <= A.n2
            A.vals[i - A.n1 + 1] = v
      else
            A.vals = push!(A.vals, zeros(typeof(A.vals[1]), i-A.n2)...)
            A.vals[end] = v
            A.n2 = i
      end   
end


function Base.:+(A::Signal, B::Signal)
      k1 = min(A.n1, B.n2)
      k2 = max(A.n2, B.n2)

      C = Signal(k1, k2, zeros(typeof(A.vals[1]), k2-k1+1))
      for i in k1:k2
            C[i] = A[i] + B[i]
      end
      return C
end

function Base.:-(A::Signal, B::Signal)
      k1 = min(A.n1, B.n2)
      k2 = max(A.n2, B.n2)

      C = Signal(k1, k2, zeros(typeof(A.vals[1]), k2-k1+1))
      for i in k1:k2
            C[i] = A[i] - B[i]
      end
      return C
end

function Base.:*(A::Signal, B::Signal)
      k1 = min(A.n1, B.n2)
      k2 = max(A.n2, B.n2)

      C = Signal(k1, k2, zeros(typeof(A.vals[1]), k2-k1+1))
      for i in k1:k2
            C[i] = A[i] * B[i]
      end
      return C
end

function convol(A::Signal, B::Signal)
      k1 = A.n1 + B.n1
      k2 = A.n2 + B.n2

      C = Signal(k1, k2, zeros(typeof(A.vals[1]), k2-k1+1))
      println((C.n1, C.n2))
      for i in k1:k2
            for j in -B.n2:-B.n1
                  C[i] += B[-j]*A[j+i]
                  println((j,j+i))
            end
            println("---")
      end
      return C
end

⋆(A::Signal, B::Signal) = convol(A::Signal, B::Signal)

# begin
#       include("signal.jl")
#       a = Signal(2,5,[1,1,1,1])
#       b = Signal(2,4,[1,1,1])
#       c = convol(a,b)
# end

# module
