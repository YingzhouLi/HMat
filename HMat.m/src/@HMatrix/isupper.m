function isu = isupper(A)

isu = false;
if A.blockType == 'H'
    for itj = 1:size(A.childHMat,2)
        for iti = itj+1:size(A.childHMat,1)
            if A.childHMat{iti,itj}.blockType ~= 'E'
                return;
            end
        end
    end
    isu = true;
end

end