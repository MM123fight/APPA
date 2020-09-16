function dataWrite(experiment,nb,nf,mi,me,n,Ae,Ai,be,bi,c)
    meta_path = [experiment, '/meta'];
    Ae_path = [experiment, '/Aeq'];
    Ai_path = [experiment, '/A'];
    be_path = [experiment, '/beq'];
    bi_path = [experiment, '/b'];
    c_path = [experiment, '/c'];
    writemeta(meta_path,nb,nf,mi,me);
    writeMat(Ai_path,Ai,mi,n);
    writeMat(Ae_path,Ae,me,n);
    writeVec(bi_path,bi);
    writeVec(be_path,be);
    writeVec(c_path,c);
end

function  writemeta(path,nb,nf,mi,me)
    name = ["nb","nf","mi","me"];
    size = [nb,nf,mi,me];
    file = fopen(path,'w');
    for row = 1:3
        fprintf(file, '%s %d\n', name(row), size(row));
    end 
    fclose(file);
    row = 4;
    file = fopen(path,'a');
    fprintf(file, '%s %d', name(row), size(row));
    fclose(file);
end

function writeMat(path,A,m,n)
    file = fopen(path,'w');
    if(m == 0)
        fprintf(file,'%d %d %.1f',m,n,0);
        fclose(file);
    else
        fprintf(file,'%d %d %.1f\n',m,n,0);
        [col,row,val] = find(A');
        s = length(row);
        if(s==1)
            fprintf(file,'%d %d %f',[row(s);col(s);val(s)]);
            fclose(file);
        else
            fprintf(file,'%d %d %f\n',[row(1:s-1)';col(1:s-1)';val(1:s-1)']);
            fclose(file);
            file = fopen(path,'a');
            fprintf(file,'%d %d %f',[row(s);col(s);val(s)]);
            fclose(file);
        end
        
    end
    
  
end

function writeVec(path,b)
    file = fopen(path,'w');
    s = length(b);
    if(s == 1)
        fprintf(file, '%f', b);
        fclose(file);
    elseif(s>1)
        fprintf(file, '%f\n', b(1:s-1));
        fclose(file);
        file = fopen(path,'a');
        fprintf(file, '%f', b(s));
        fclose(file);
    end 
end