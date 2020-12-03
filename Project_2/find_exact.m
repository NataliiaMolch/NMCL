% Function evaluates the exact solution at a given time

function Uexact = find_exact(pIC,vIC,S,Sinv,Lambda,xc,time)

VSelect = @(V,IDX) V(IDX,:);
V0 =@(x) Sinv*[pIC(x);vIC(x)];
VT_1 =@(x) VSelect(V0(x),1);
VT_2 =@(x) VSelect(V0(x),2);

Uexact = S*[VT_1(xc-Lambda(1,1)*time);
            VT_2(xc-Lambda(2,2)*time)];
        
        