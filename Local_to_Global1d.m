function gp=Local_to_Global1d(lp,gv)

gp=zeros(size(gv,1),size(lp,1));

for i=1:length(lp)
    gp(:,i) = (gv(:,1)+gv(:,2))/2.0+lp(i)*(gv(:,2)-gv(:,1))/2.0;
end

end