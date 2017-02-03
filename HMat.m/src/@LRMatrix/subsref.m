function B = subsref(A,S)

if numel(S) > 1
    B = builtin('subsref',A,S);
    return;
end

switch S(1).type
    case '()'
        rows = S(1).subs{1};
        cols = S(1).subs{2};
        B = LRMatrix(A.UMat(rows,:),A.VMat(cols,:),A.EPS,A.MAXRANK);
    case '{}'
        rows = S(1).subs{1};
        cols = S(1).subs{2};
        if length(rows) == 1 && length(cols) == 1
            B = A.UMat(rows,:)*A.VMat(cols,:)';
        else
            B = builtin('subsref',A,S);
        end
    otherwise
        B = builtin('subsref',A,S);
end

end