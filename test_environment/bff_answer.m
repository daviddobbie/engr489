%% BFF Answer
% David Dobbie
% Victoria University of Wellington

function [answer] = bff_answer(f_answer, g_bfv)
    
    g_poro = ones(length(f_answer),1);
    bfv = g_bfv' * f_answer;
    poro = g_poro' * f_answer;
    
    answer = bfv/poro;
end

