function pale=assembly_mix_with_white(color,amount)
%ASSEMBLY_MIX_WITH_WHITE Blend an RGB color toward white.
pale=color+(1-color).*amount;
end
