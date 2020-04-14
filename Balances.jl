using DelimitedFiles

A = readdlm("C:\\Users\\user\\Desktop\\Cornell Academics\\Spring 2020\\CHEME 7770\\PS3\\Atom_Matrix.dat")
display(A)
size(A)

S = readdlm("C:\\Users\\user\\Desktop\\Cornell Academics\\Spring 2020\\CHEME 7770\\PS3\\S_Matrix.dat")
display(S)
size(S)

res = transpose(A)*S
display(res)
(nrow,ncol) = size(res)

s = zeros(nrow)
for i in 1:nrow
    for j in 1:ncol
        s[i] = s[i] + res[i,j]
    end
end

display(s)
