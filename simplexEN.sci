function simplex(A,b,rel,c,typ)
    
    tic();
    //A : coefficients matrix
    //n : column vector of constraints
    //rel : column vector precising the type of constraint (-1: <= / 0: == / 1: >=)
    //c : row costs vector 
    //typ : type of objective : "min" or "max"
    
    nbVar = length(c)
    
    // z(opt) is initialized at 0
    z=0
    
    // Computing and displaying the problem under its standard form
    display("STANDARD FORM")
    [A,b,c,base]=standardForm(A,b,rel,c,typ);
    disp(A); disp(b); disp(c); disp(z); disp(base)
    
    // Le vecteur base est un vecteur ligne. Il contient l'information des vecteurs qui sont en base ainsi que leur position (ligne). Ainsi, si le vecteur base est (2,3,5), la première ligne de A a pour vecteur en base le vecteur 2, la seconde ligne de A a pour vecteur en base le vecteur 3...
    
    if isBaseValid(A,base)==1 then
        
        // if the problem is already in a standard form, the base is feasible, simplex is directly applied.
        display("NORMAL ITERATION")
        
        [A,b,z,c,null,null,base]=iterSimplex(A,b,z,c,null,null, base)
        
        display("FINAL RESULT")
        displayTab(A,b,c,z,base,nbVar)
        
        listsol = degeneracyv2(A,b,c,z,base,nbVar, list(getSolution(b, base, nbVar)))
        displaydegeneracy(listsol)
    else
        
        // if turning the problem into a standard form doesn't return a feasible base, it is processed through the phase 1
        //...then simplex is performed
        [ok,A,b,c,base,z]=phaseOne(A,b,c,base)
        if ok==1 then
            display("ITERATION CLASSIQUE")
            
            [A,b,z,c,null,null,base]=iterSimplex(A,b,z,c,null,null, base)
            
            display("FINAL RESULT")
            displayTab(A,b,c,z,base,nbVar)
            
            listsol = degeneracyv2(A,b,c,z,base,nbVar, list(getSolution(b, base, nbVar)))
            displaydegeneracy(listsol)
        end
    end
    
    // displaying runtime
    display("Runtime (in seconds) : ")
    disp(toc())
    
endfunction





function display(text)
    disp("-----------------------------------------------------------------------")
    disp(text)
    disp("-----------------------------------------------------------------------")
endfunction





function [A,b,c,base]=standardForm(A, b, rel, c, typ)
    
    prevNbColA = size(A,2)
    
    // Consistency of matrix dimensions is checked on rows
    if ~isequal([size(A,1), size(b,1), size(rel,1)]/size(A,1), [1,1,1]) | (size(c,1)~=1) then
        disp("Number of rows is inconsistent on some matrices")
        return null;
    end
    
    // Consistency of matrix dimensions is checked on columns
    if ~isequal([size(A,2), size(c,2)]/size(A,2), [1,1]) | ~isequal([size(b,2), size(rel,2)], [1,1]) then
        disp("Number of columns is inconsistent on some matrices")
        return null;
    end
    
    // stack variables are generated
    stackVar=eye(size(A,1), size(A,1))
    for i=1:size(rel,1)
        stackVar(i,i)=stackVar(i,i)*(-1)*rel(i)
    end
    
    // The base(now feasible) that was just created is inilialized
    base = prevNbColA+1:prevNbColA+size(stackVar,2)


    // basic variables are added to A (if non null). 
    // if the stack variable shouldn't be integrated, the base is corrected   
    for i=1:size(stackVar,2)
        if ~isequal(stackVar(:,i), zeros(size(stackVar,1),1)) then
            A(:, size(A,2)+1)=stackVar(:, i)
            c(:, size(c,2)+1)=0
        else
            base(i)=-1
        end
    end
    
    // if the cost function is to minimize, c is multiplied by -1
    if typ=="min" then
        c=c*(-1)
    end
    
    // the row equations are arranged so that all b elements are positive
    for i=1:size(b,1)
        if b(i)<0 then
            b(i) = b(i)*(-1)
            A(i,:) = A(i,:)*(-1)
        end
    end
    
endfunction



function isOk=isBaseValid(A,base)
    isOk=1
    for i=base
        if i==-1 | A(find(base==i,1),i)==-1 then
            isOk=0
            return isOk
        end
    end
    return isOk
endfunction



function [ok,A,b,c,base,z]=phaseOne(A,b,c,base)
    
    ok=1
    
    display("PHASE 1")
    // AP1 and AC1 are used for phase 1
    Ap1 = A
    cp1 = zeros(1,length(c))
    zp1 = 0
    
    // Used to merge A and auxiliary variables 
    I=eye(size(A,1),size(A,1))
    
    // for each element of the base eventual sources of error are checked
    for i=1:length(base)
        if base(i)==-1 | Ap1(i,base(i))==-1 then
            Ap1(:,$+1)=I(:,i)
            cp1(1,$+1)=0
            cp1(1,1:length(c))=cp1(1,1:length(c))+Ap1(i,1:length(c))
            zp1 = zp1 + b(i)
            base(i)=size(Ap1, 2)
        end
        
    end
    
    cp1bis = cp1
    cp1bis(1:length(c))=c
    
    disp(Ap1)
    disp(cp1)
    disp(zp1)
    
    // Simplex iteration is applied to the phase 1
    //the initial z and cost vector c are also passed to the function so that they are kept updated
    [Ap1,b,zp1,cp1,z,cp1bis,base]=iterSimplex(Ap1,b,zp1,cp1,z,cp1bis, base)
    
    if zp1~=0 then
        disp("Phase 1 cannot delete auxiliary variables : no finite solution to the phase 1 problem")
        ok=-1
    end
    
    // Auxiliary variables are deleted from A and c before the end of phase 1
    // A,b,c,z are then ready for simplex iterations (if the actual base is feasible)
    A=Ap1(:,1:size(A,2))
    c=cp1bis(1:length(c))
    
endfunction



function [A,b,z1,c1,z2,c2,base]=iterSimplex(A,b,z1,c1,z2,c2,base)
    
    // Iterations are performed while there is a s>0
    while chooseS(c1)~=(-1)
        // Choosing r
        in = chooseS(c1)
        // Choosing s
        out = chooseR(A,b,in)
        // If r not found (Ar<=0), no finite solution
        if out==(-1) then
            disp("Error, pas de solution finie")
            return
        end
        // Exchanging in and out variables in the base
        base(out)=in
        // pivot is performed
        [A,b,z1,c1,z2,c2]=pivot(A,b,z1,c1,z2,c2,in,out,base)
        
        disp("_________________")
        disp("Next iteration (A,b,c,z)")
        disp(A); disp(b); disp(c1); disp(z1)
    end
    
endfunction




// Choosing the set column
function [in]=chooseS(c, base)
    maxi=max(c)
    if maxi<=0 then
        in=-1
    else
        in = find(c==maxi,1)
    end
endfunction



// Choosing the reset line
function [out]=chooseR(A,b,in)
    out=-1
    mini=-1
    Ar=A(:,in)
    for i=1:length(Ar)
        if (Ar(i)>0 & ((b(i)/Ar(i))<mini | mini==-1)) then
            mini=b(i)/Ar(i)
            out=i
        end
    end
endfunction




// Pivot operation function
function [A,b,z1,c1,z2,c2]=pivot(A,b,z1,c1,z2,c2,in,out,base)
    Ai=A(:,base)
    
    z1=z1-(c1(in)/A(out,in))*b(out)
    if (z2~=null) then
        z2=z2-(c2(in)/A(out,in))*b(out)
    end
    c1=c1-(c1(in)/A(out,in))*A(out,:)
    if (c2~=null) then
        c2=c2-(c2(in)/A(out,in))*A(out,:)
    end
    
    A=Ai\A
    b=Ai\b
    
endfunction


function xfin=getSolution(b,base,nbVar)
    xfin = zeros(nbVar,1);
    for i=1:length(xfin)
        if ~isempty(find(base==i)) then
            xfin(i)=b(find(base==i))
        end
    end
endfunction




// display des états respectifs de A, b, c, z à un instant donné
function displayTab(A,b,c,z,base,nbVar)
    disp("Valeur de A ...");
    disp(A); 
	disp("... de b ...");
	disp(b); 
	disp("... de c ...");
	disp(c); 
	disp("... de z ...");
	disp(z); 
	disp("Vecteur des variables de decision ");
	disp(getSolution(b,base,nbVar));
endfunction





function listsol=degeneracyv2(A,b,c,z,base,nbVar, listsol)
    // On cherche une eventuelle dégénérescence de seconde espèce
    in=1
    totest=list()
    for in=1:length(c)
        if isempty(find(base==in)) & c(in)==0 then
            totest($+1)=in
        end
        in=in+1
    end
    
    // On a alors stocké dans la liste 'totest' toutes les positions de 0 dans le vecteur coût (hors base) qui sont probablement sinonymes de nouvelle solution valide
        
    oldbase=base
    oldb=b
    oldz=z
    oldc=c
    oldA=A
       
    // Pour chaque position de 0 repéré dans la liste 'totest', on effectue une nouvelle itération du simplexe qui nous donne possiblement une solution. Si solution il y a, on va tester si cette solution est nouvelle ou non. Si elle l'est, on va aller en profondeur voir si on ne trouve pas de nouvelle solution à partir de cette dernière.
    for idx=1:length(totest)
        base=oldbase
        b=oldb
        z=oldz
        c=oldc
        A=oldA
        
        in = totest(idx)
        
        out = chooseR(A,b,in)
        base(out)=in
        [A,b,z,c]=pivot(A,b,z,c,null,null,in,out,base)
        
        newsol=getSolution(b,base,nbVar)

        if myFind(listsol,newsol)==0 then
            listsol($+1)=newsol
            listsol=degeneracyv2(A,b,c,z,base,nbVar, listsol);
        end
    end
    
endfunction



// Fonction permettant de voir si une solution passée en argument est déjà solution (ou non) d'une liste de solutions également passée en argument
function isIn=myFind(listsol, sol)
    isIn = 0
    for i=1:length(listsol)
        if listsol(i)==sol then
            isIn = 1
        end
    end
endfunction



// Fonction d'display de la dégénérescence (pour une liste de solutions passée en argument)
function displaydegeneracy(listsol)
    
    if length(listsol)>1 then
        display("DEGENERESCENCE DE SECONDE ESPECE")
        disp("Les solutions sont les combinaisons convexes de...")
        disp(listsol(1))
        for idx=2:length(listsol)
            disp("...et...")
            disp(listsol(idx))
        end
    end
    
endfunction
