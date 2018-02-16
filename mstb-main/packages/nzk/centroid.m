%Calculate the centroid of two vectors aT, bT
%result on aT
function cT = centroid(aT,bT)
cT = dot(aT , bT)/sum(bT);
end