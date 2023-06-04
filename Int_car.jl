using Plots
using LinearAlgebra
#plotlyjs()
#gr(dpi=450)
size=(dpi=400)
#Code to convert from internal coordintes to cartesian coordinates
##
function Int_Cart(R,θ,α,β)
    rH = 0.74#H2 bond distance
    rC = [1.26193, 1.34953, 1.23481, 1.36188, 1.17254]#C-C and C-N bond distances
    Coef = 4:-1:1 #Coeff to caclulate the bond deistances
    CoM_CN = (12*sum(Coef.*rC[1:4]) + 14*sum(rC))/(4*12+14) #centeer of mass of C5N group
    C_cord = [sum(rC[1:i]) for i in 1:5] 
    pushfirst!(C_cord,0) #The first carbon atom position
    C_cord = C_cord .- CoM_CN # We displace the origin so that the CoM(C5N) is at the orgin
    Cord_mat = zeros(8,3)
    Cord_mat[1:6,1] = C_cord
    H_CoM = [R*cosd(θ),R*sind(θ),0]
    r = rH/2 #distance from CoM(H2) to H
    R_cap = [cosd(θ),sind(θ),0]
    R_per_Cap = [sind(θ),-cosd(θ),0]
    r_h1 = (r*cosd(α).*R_cap + r*sind(α)*cosd(β).*R_per_Cap + [0,0,r*sind(α)*sind(β)]) + H_CoM
    r_h2 =  H_CoM - (r*cosd(α).*R_cap + r*sind(α)*cosd(β).*R_per_Cap + [0,0,r*sind(α)*sind(β)])  #r*cosd(180+α).*R_cap + r*sind(180+α)*cosd(180+β).*R_per_Cap + [0,0,r*sind(180+α)*sind(180+β)]
    R_H = vcat(r_h1',r_h2')
    Cord_mat[7:8,:] = R_H
    return Cord_mat
end
C = ["black","black","black","black","black","blue","grey","grey"]
##
Cord=Int_Cart(4,90,90,90)
scatter(Cord[:,1],Cord[:,2],Cord[:,3],color=C, xlimits=(-5,5), ylimits=(-5,5), zlimits=(-2,2))
##
R = LinRange(1,4,100)
Th = LinRange(60,120,100)
al = LinRange(0,180,300)
bt = LinRange(0,180,300)
Rr = zeros(200)
anim = @animate for i in 1:300
    j = Int(floor(i/3))
    if j==0
        j=1
    end
    Cord=Int_Cart(2.5,90,90,bt[i])
    #ϵ=(sum((Cord[7,:] - Cord[8,:]).^2))^0.5
    #Rr[i] = ϵ
    scatter(Cord[1:6,1],Cord[1:6,2],Cord[1:6,3],color=C[1:6], xlimits=(-4,4), ylimits=(-3,3), zlimits=(-4,4),markersize=10)
    scatter!(Cord[7:8,1],Cord[7:8,2],Cord[7:8,3],color=C[7:8],markersize=8,xaxis=false,yaxis=false,zaxis=false)
    plot!(Cord[7:8,1],Cord[7:8,2],Cord[7:8,3],color=C[7:8],linewidth=4,camera = (60, 15),legend=false)
end
##
gif(anim, "C5N_H2.gif", fps = 30)
##
#plot([1.26193,0],[],color=["grey"])
###
#r=(Cord[7,:] - Cord[8,:])
###
#Cord=Int_Cart(2.5,50,30,60)
#scatter(Cord[1:6,1],Cord[1:6,2],Cord[1:6,3],color=C[1:6], xlimits=(-4,4), ylimits=(-3,3), zlimits=(-4,4),markersize=10)
#scatter!(Cord[7:8,1],Cord[7:8,2],Cord[7:8,3],color=C[7:8],markersize=8,xaxis=false,yaxis=false,zaxis=false)
#plot!(Cord[7:8,1],Cord[7:8,2],Cord[7:8,3],color=C[7:8],linewidth=4,camera = (60, 10),legend=false)
#savefig("Dummy.png")
##
#Vel = [] #Electronic repulsion vector
#q = [6,6,6,6,6,7,1,1]#charge vector
#for i in 1:7
#    for j in i+1:8
#        vel = q[i]*q[j]/norm(Cord[i]-Cord[j])
#        append!(Vel,vel)
#    end
#end
#Vel = sum(Vel)