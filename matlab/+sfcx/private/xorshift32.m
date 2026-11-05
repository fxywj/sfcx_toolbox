function out = xorshift32(seed)
%XORSHIFT32 Deterministic xorshift32 PRNG for cross-language reproducibility.
%   Returns a function handle that generates the next random uint32.
    x = uint32(seed);
    out = @next;
    function r = next()
        x = bitxor(x, bitshift(x, 13));
        x = bitxor(x, bitshift(x, -17));
        x = bitxor(x, bitshift(x, 5));
        r = x;
    end
end
