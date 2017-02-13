function isl = islower(A)

isl = false;
if A.blockType == 'H'
    for iti = 1:size(A.childHMat,1)
        for itj = iti+1:size(A.childHMat,2)
            if A.childHMat{iti,itj}.blockType ~= 'E'
                return;
            end
        end
    end
    isl = true;
end

end