function [nb,nf,mi,me,m,n,Ae,Ai,be,bi,c] = dataRead(experiment)
    meta_path = [experiment, '/meta'];
    Ae_path = [experiment, '/Aeq'];
    Ai_path = [experiment, '/A'];
    be_path = [experiment, '/beq'];
    bi_path = [experiment, '/b'];
    c_path = [experiment, '/c'];
    
    %read data
    [nb,nf,mi,me] = readmeta(meta_path);
    m = mi + me;
    n = nb + nf;
    Ae = readMat(Ae_path,me,n);
    Ai = readMat(Ai_path,mi,n);
    be = readVec(be_path,me);
    bi = readVec(bi_path,mi);
    c = readVec(c_path,n);
end

function [nb,nf,mi,me] = readmeta(path)
    file = fopen(path);
    size = textscan(file,'%s %n');
    nb = size{2}(1);
    nf = size{2}(2);
    mi = size{2}(3);
    me = size{2}(4);
    fclose(file);
end

function [A] = readMat(path,m,n)
    file = fopen(path);
    A_data = textscan(file,'%n%n%f');
    if(m ~= A_data{1}(1))
        error(['row size is not accordant in ', path]);
    end
    if(n ~= A_data{2}(1))
        error(['col size is not accordant in ', path]);
    end
    A_data{1}(1) = [];
    A_data{2}(1) = [];
    A_data{3}(1) = [];
    A = sparse(A_data{1},A_data{2},A_data{3},m,n);
    fclose(file);
end

function [b] = readVec(path,n)
    file = fopen(path);
    b_data = textscan(file,'%f');
    b = b_data{1};
    if(length(b) ~= n)
        error(['vector size is not accordant in ', path]);
    end
    fclose(file);
end


