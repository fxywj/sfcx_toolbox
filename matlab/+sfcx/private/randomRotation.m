function R = randomRotation(rng)
%RANDOMROTATION Uniform random rotation matrix from a unit quaternion.
q = randn(rng, 4, 1);
q = q / norm(q);
w = q(1); x = q(2); y = q(3); z = q(4);
R = [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w); ...
     2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w); ...
     2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)];
end
