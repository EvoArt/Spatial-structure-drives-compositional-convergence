# Alpha
function shannon(a)
    a .= [isnan(x) ? 0.0 : x for x in a]
    a = a[a .>0.0]
    A = sum(a)
    H = -sum([(a[i]/A) * log(a[i]/A) for i in 1:length(a)])
    return isnan(H) ? 0 : H
end

function pielou(a) 
    hmax = shannon(a) / log(sum(a .> 0))
    return isnan(hmax) ? 0 : hmax
end

# Beta
function βJAC(X)
    x = Bool.(X)
    n_sites = size(x)[1]
    b(x₁,x₂) = sum(x₁ .& .!x₂)
    num = 0
    denom = sum(x[1,:])-(sum(sum(x,dims = 1) .>0)) 
    for i in 2:n_sites
        denom +=(sum(x[i,:]))

        for j in 1:i
            num  += b(x[i,:],x[j,:]) 
            num += b(x[j,:],x[i,:])
        end
    end
    num / (num +denom)
end

function βJTU(X)
    x = Bool.(X)
    n_sites = size(x)[1]
    b(x₁,x₂) = sum(x₁ .& .!x₂)
    num = 0
    denom = sum(x[1,:])-(sum(sum(x,dims = 1) .>0)) 
    for i in 2:n_sites
        denom +=(sum(x[i,:]))
        for j in 1:i
            num  += 2min(b(x[i,:],x[j,:]), b(x[j,:],x[i,:]))
        end
    end
    num / (num +denom)
end
βJNE(X) = βJAC(X) - βJTU(X)
βpart(X,sigfig = 2) = (βJAC = round(βJAC(X),digits = sigfig),βJTU = round(βJTU(X),digits = sigfig), βJNE = round(βJNE(X),digits = sigfig),Relative_turnover = round(βJTU(X)/ βJAC(X),digits = sigfig))

