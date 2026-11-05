function z = fisherZ(r)
%FISHERZ Fisher z-transform, clipped to [-0.9999, 0.9999].
z = atanh(max(min(r, 0.9999), -0.9999));
end
