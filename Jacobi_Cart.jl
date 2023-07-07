using Plots
using LinearAlgebra
#plotlyjs()
#gr(dpi=450)
gr(size=(500,500))
#Code to convert from internal coordintes to cartesian coordinates
#R distance between Center of Mass (CoM) of 2 fragments
#r1, r2 internuclear distances of 2 fragments
# θ,α,β are the angles difined by standard Jacobi.
#θ the angle between First fragment and CoM of second fragment.
#α is the angle between the second fragment and the CoM of the first fragment
##
function Int_Cart_C5N(R,θ,α,β;rH = 0.74)#For C5N + H2 system
    #rH is the H2 bond distance
    #All angles in degree
    rC = [1.26193, 1.34953, 1.23481, 1.36188, 1.17254]#C-C and C-N bond distances
    Coef = 4:-1:1 #Coeff to caclulate the bond deistances
    CoM_CN = (12*sum(Coef.*rC[1:4]) + 14*sum(rC))/(5*12+14) #centeer of mass of C5N group
    C_cord = [sum(rC[1:i]) for i in 1:5] 
    pushfirst!(C_cord,0) #The first carbon atom position
    C_cord = C_cord .- CoM_CN # We displace the origin so that the CoM(C5N) is at the orgin
    Cord_mat = zeros(8,3)#The coordinate matrix
    Cord_mat[1:6,1] = C_cord#Copying the xcomponents to the coordinate matrix
    H_CoM = [R*cosd(θ),R*sind(θ),0]
    r = rH/2 #distance from CoM(H2) to H
    R_cap = [cosd(θ),sind(θ),0]#The vector from the CoM of the first fragment to the second fragment
    R_per_Cap = [sind(θ),-cosd(θ),0]
    r_h1 = (r*cosd(α).*R_cap + r*sind(α)*cosd(β).*R_per_Cap + [0,0,r*sind(α)*sind(β)]) + H_CoM
    r_h2 =  H_CoM - (r*cosd(α).*R_cap + r*sind(α)*cosd(β).*R_per_Cap + [0,0,r*sind(α)*sind(β)])  #r*cosd(180+α).*R_cap + r*sind(180+α)*cosd(180+β).*R_per_Cap + [0,0,r*sind(180+α)*sind(180+β)]
    R_H = vcat(r_h1',r_h2')
    Cord_mat[7:8,:] = R_H
    return Cord_mat
end
function Jacobi_Cart(R,θ,α,β;r1=1,r2=1, m1 = [1,1], m2 = [2,2])# For general Jacobi with 4 atoms
    #All angles in degree
    Cord_mat = zeros(4,3) #the coordinate matrix; 4 atoms with x,y, z coordinates
    #first place the first atom of the first fragment at origin and the second at the r1
    Frag1 = [0,r1] #This is just the x coordinates
    # center of mass of first fragment in at origin
    CoM_Frag1 = r1*m1[2]/sum(m1) # the CoM equation
    Frag1 = Frag1 .- CoM_Frag1 #Placing the CoM at origin
    Cord_mat[1:2,1] = Frag1
    CoM_Frag2 = [R*cosd(θ),R*sind(θ),0]#second fragments CoM is at the xy plane
    rm1 = r2*m2[2]/(sum(m2))#Distance from CoM of frag2 and first mass of frag2
    rm2 = -r2*m2[2]/(sum(m2))#Negative sign because it is in the oppposite direction t0 mass 1
    R_cap = [cosd(θ),sind(θ),0]#The vector from the CoM of the first fragment to the second fragment
    R_per_Cap = [sind(θ),-cosd(θ),0]
    rm1_vec = (rm1*cosd(α).*R_cap + rm1*sind(α)*cosd(β).*R_per_Cap + [0,0,rm1*sind(α)*sind(β)]) +  CoM_Frag2 # The position vector of first mass of fragment 2
    rm2_vec = (rm2*cosd(α).*R_cap + rm2*sind(α)*cosd(β).*R_per_Cap + [0,0,rm2*sind(α)*sind(β)]) +  CoM_Frag2 # The position vector of second mass of fragment 2
    Frag2 = vcat(rm1_vec', rm2_vec')#concantenating the vectors
    Cord_mat[3:4,:] = Frag2#Adding to the coordinate matrix
    return Cord_mat
end
C = ["black","black","black","black","black","blue","grey","grey"]
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
    Cord=Int_Cart_C5N(2.5,90,90,bt[i])
    #ϵ=(sum((Cord[7,:] - Cord[8,:]).^2))^0.5
    #Rr[i] = ϵ
    scatter(Cord[1:6,1],Cord[1:6,2],Cord[1:6,3],color=C[1:6], xlimits=(-4,4), ylimits=(-3,3), zlimits=(-4,4),markersize=10)
    scatter!(Cord[7:8,1],Cord[7:8,2],Cord[7:8,3],color=C[7:8],markersize=8,xaxis=false,yaxis=false,zaxis=false)
    plot!(Cord[7:8,1],Cord[7:8,2],Cord[7:8,3],color=C[7:8],linewidth=4,camera = (60, 15),legend=false)
end
##
gif(anim, "C5N_H2.gif", fps = 30)
##
r = LinRange(0.5,1.5,50)
R = LinRange(1,4,100)
Th = LinRange(60,120,400)
al = LinRange(0,180,400)
bt = LinRange(0,180,400)
Rr = zeros(200)
n = vcat(1:50,50:-1:1)
n = vcat(n,n)
n = vcat(n,n)
anim = @animate for i in 1:400
    r1 = r[n[i]]
    r2 = r[end-n[i]+1]
    Cord=Jacobi_Cart(2.0,Th[i],al[i],90;r1=r1,r2=r2)
    #ϵ=(sum((Cord[7,:] - Cord[8,:]).^2))^0.5
    #Rr[i] = ϵ
    scatter(Cord[1:4,1],Cord[1:4,2],Cord[1:4,3],color="black", xlimits=(-4,4), ylimits=(0,4), zlimits=(-4,4),markersize=10)
    plot!(Cord[1:2,1],Cord[1:2,2],Cord[1:2,3],color="grey",linewidth=3,camera = (60, 15),legend=false)
    plot!(Cord[3:4,1],Cord[3:4,2],Cord[3:4,3],color="grey",linewidth=3)
end
##
gif(anim, "H2_H2.gif", fps = 30)
