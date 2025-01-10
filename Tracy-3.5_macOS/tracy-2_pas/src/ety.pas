  module ety(input, output);

  %include 'mathlib.def'
  %include 'ety.def'

  { Interface }

  [global] procedure ETY(n, low, igh : integer; var a : matrix; var ort : vector); forward;

  [global] procedure ETYT(n, low, igh : integer; var a : matrix; var ort : vector;
		 var z : matrix); forward;

  [global] procedure ety2(n, low, igh : integer; var h : matrix; var wr, wi : vector;
		 var z : matrix; var ierr : integer); forward;

  { Implementation }

    procedure ETY{n, low, igh : integer; var a : matrix; var ort : vector};

      var
           i,j,m,ii,jj,la,mp,kp1:integer;
           f,g,h,scale:double;

    function dsign(x,y:double):double;
      begin
        x:=abs(x);
        if y>=0 then dsign:=x
                else dsign:=-x;
      end;

{    this subroutine is a translation of the algol procedure orthes,
     num. math. 12, 349-368(1968] by Martin and Wilkinson.
     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971].

     given a real general matrix, this subroutine
     reduces a submatrix situated in rows and columns
     low through igh to upper hessenberg form by
     orthogonal similarity transformations.

     on input-

        n is the order of the matrix,

        low and igh are integers determined by the balancing
          subroutine  balanc.  ifbalanc  has not been used,
          set low:=1, igh:=n,

        a contains the input matrix.

     on output-

        a contains the hessenberg matrix.  information about
          the orthogonal transformations used in the reduction
          is stored in the remaining triangle under the
          hessenberg matrix,

        ort contains further information about the transformations.
          only elements low through igh are used.

     fortran routine by b. s. garbow
     modified by filippo neri. }

begin
  la:= igh-1;
  kp1:= low+1;
     { if la < kp1 goto 200}
  if la>= kp1 then
    begin
      for {180} m:= kp1 to la do
      begin
         h:= 0.0;
         ort[m]:= 0.0;
         scale:= 0.0;

{     ********** scale column (algol tol then not needed] **********}

         for {90} i:= m to igh do scale:= scale+Abs(a[i,m-1]);

        { if scale = 0.0 goto 180 }

         if scale<>0 then
         begin
           mp:= m+igh;
{     ********** for i:=igh step -1 until m for -- **********}

           for {100} ii:= m to igh do
           begin
              i:= mp-ii;
              ort[i]:= a[i,m-1] / scale;
              h:= h+ort[i]*ort[i];
    {100}  end;

           g:= -dsign(sqrt(h),ort[m]);
           h:= h-ort[m]*g;
           ort[m]:= ort[m]-g;
{     ********** form (i-(u*ut]/h]*a **********}
         for {130} j:= m to n do
         begin
            f:= 0.0;
{     ********** for i:=igh step -1 until m for -- **********}
            for {110} ii:= m to igh do
            begin
               i:= mp-ii;
               f:= f+ort[i]*a[i,j];
  {110}     end;

            f:= f / h;

            for {120} i:= m to  igh do
    {120}       a[i,j]:= a[i,j]-f*ort[i];

  {130}  end;

{     ********** form (i-(u*ut]/h]*a*(i-(u*ut]/h] **********}
         for {160} i:= 1 to igh do
         begin
            f:= 0.0;
{     ********** for j:=igh step -1 until m for -- **********}
            for {140} jj:= m to igh do
            begin
               j:= mp-jj;
               f:= f+ort[j]*a[i,j];
  {140}     end;

            f:= f / h;

            for j:= m to igh do a[i,j]:= a[i,j]-f*ort[j];

  {160}  end;

         ort[m]:= scale*ort[m];
         a[m,m-1]:= scale*g;
      end;   
  {180} end;
   end;
 end;

  procedure ETYT{n, low, igh : integer; var a : matrix; var ort : vector;
	         var z : matrix};

  var i,j,kl,mm,mp,mp1	: integer;
      g			: double;

{    this subroutine is a translation of the algol procedure ortrans,
     num. math. 16, 181-204[1970] by peters and wilkinson.
     handbook for auto. comp., vol.ii-linear algebra, 372-395[1971].

     this subroutine accumulates the orthogonal similarity
     transformations used in the reduction of a real general
     matrix to upper hessenberg form by  ety.

     on input-

        n is the order of the matrix,

        low and igh are integers determined by the balancing
          subroutine  balanc.  if  balanc  has not been used,
          set low:=1, igh:=n,

        a contains information about the orthogonal trans-
          formations used in the reduction by  orthes
          in its strict lower triangle,

          ort contains further information about the trans-
          formations used in the reduction by  ety.
          only elements low through igh are used.

     on output-

        z contains the transformation matrix produced in the
          reduction by  ety,

        ort has been altered.

     fortran routine by b. s. garbow.
     modified by f. neri.


     ********** initialize z to identity matrix **********}

begin
     for i := 1 to n do
      begin
         for j := 1 to n do
             z[i,j] := 0.0;
             z[i,i] := 1.0;
      end;

    kl := igh - low - 1;
    if kl >= 1 then
    begin
      for mm := 1 to kl do
      begin
        mp := igh - mm;
        if a[mp,mp-1] <> 0.0 then
        begin
          mp1 := mp + 1;
          for i := mp1 to igh do
                 ort[i] := a[i,mp-1];
          for  j := mp to igh do
          begin
            g := 0.0;
            for i := mp to igh do
                g := g + ort[i] * z[i,j];
{     ********** divisor below is negative of h formed in orthes.
                 double division avoids possible underflow **********}
            g := (g / ort[mp]) / a[mp,mp-1];
            for i := mp to igh do
               z[i,j] := z[i,j] + g * ort[i];
          end;
       end;
      end;
    end;
  end;

  procedure ety2{n, low, igh : integer; var h : matrix; var wr, wi : vector;
		 var z : matrix; var ierr : integer};

{    this subroutine is a translation of the algol procedure hqr2,
     num. math. 16, 181-204[1970] by peters and wilkinson.
     handbookfor{auto. comp., vol.ii-linear algebra, 372-395[1971].

     this subroutine finds the eigenvalues and eigenvectors
     of a real upper hessenberg matrix by the qr method.  the
     eigenvectors of a real general matrix can also be found
     ifelmhes  and  eltran  or  orthes  and  ortran  have
     been used to reduce this general matrix to hessenberg form
     and to accumulate the similarity transformations.

     on input-

        n is the order of the matrix,

        low and igh are integers determined by the balancing
          subroutine  balanc.  ifbalanc  has not been used,
          set low:=1, igh:=n,

        h contains the upper hessenberg matrix,

        z contains the transformation matrix produced by  eltran
          after the reduction by  elmhes, or by  ortran  after the
          reduction by  orthes, ifperformed.  ifthe eigenvectors
          of the hessenberg matrix are desired, z must contain the
          identity matrix.

     on output-

        h has been destroyed,

        wr and wi contain the real and imaginary parts,
          respectively, of the eigenvalues.  the eigenvalues
          are unordered except that complex conjugate pairs
          of values appear consecutively with the eigenvalue
          having the positive imaginary part first.  ifan
          error exit is made, the eigenvalues should be correct
         for{indices ierr+1,...,n,

        z contains the real and imaginary parts of the eigenvectors.
          ifthe i-th eigenvalue is real, the i-th column of z
          contains its eigenvector.  ifthe i-th eigenvalue is complex
          with positive imaginary part, the i-th and [i+1]-th
          columns of z contain the real and imaginary parts of its
          eigenvector.  the eigenvectors are unnormalized.  ifan
          error exit is made, none of the eigenvectors has been found,

        ierr is set to
          zero      for{normal return,
          j          ifthe j-th eigenvalue has not been
                     determined after 30 iterations.

     arithmetic is double precision. complex division
     is simulated by routin etdiv.

     fortran routine by b. s. garbow.
     modified by f. neri.


     ********** machep is a machine dependent parameter specifying
                the relative precision of floating point arithmetic.

                **********}

      label 60,70,100,150,260,270,280,320,330,340,1000;

      const	machep=1e-20;

      var i,j,k,l,m,en,ii,jj,ll,mm,na,nn:integer;
          its,mp2,enm2:integer;
          p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm:double;
          notlas:boolean;
          z3r,z3i:double;

    function dsign(x,y:double):double;
      begin
        x:=abs(x);
        if y>=0 then dsign:=x
                else dsign:=-x;
      end;


  function min0(i,j:integer):integer;
    begin
      if i>=j then min0:=j 
              else min0:=i
    end;


  procedure etdiv(var a,b:double; c,d,e,f:double);
{---------------------------------------------------
   computes the complex division
     a + ib := (c + id)/(e + if)
  very slow, but tries to be as accurate as
  possible by changing the order of the
  operations, so to avoid under(over)flow
  problems.
  Written by F. Neri Feb. 12 1986
 --------------------------------------------------}

   var s,t,cc,dd,ee,ff,temp:double;
       flip:integer;
   begin

      flip := 0;
      cc := c;
      dd := d;
      ee := e;
      ff := f;
      if  Abs(f)>=Abs(e)  then
      begin
        ee := f;
        ff := e;
        cc := d;
        dd := c;
        flip := 1;
      end;
      s := 1.0/ee;
      t := 1.0/(ee+ ff*(ff*s));
      if  Abs(ff) >= Abs(s)  then
      begin
        temp := ff;
        ff := s;
        s := temp;
      end;

      if  Abs(dd) >= Abs(s)  then
        a := t*(cc + s*(dd*ff))
      else if  Abs(dd) >= Abs(ff)  then
        a := t*(cc + dd*(s*ff))
      else
        a := t*(cc + ff*(s*dd));
 
      if  Abs(cc) >= Abs(s) then
        b := t*(dd - s*(cc*ff))
      else if  Abs(cc) >= Abs(ff) then
        b := t*(dd - cc*(s*ff))
      else
        b := t*(dd - ff*(s*cc));

      if flip<>0  then b := -b;

   end;

begin
      ierr:=0;
      norm:=0.0;
      k:=1;
{     ********** store roots isolated by balanc
                 and compute matrix norm ********** }
     for{50} i:=1 to n do
     begin
        for{40} j:=k to n do
   {40}    norm:=norm+ABS(h[i,j]);
           k:=i;
      {   if (i >= low) and (i <= igh) goto 50}
           if not((i >= low) and (i <= igh)) then
           begin
             wr[i]:=h[i,i];
             wi[i]:=0.0;
           end;
{50} end;

      en:=igh;
      t:=0.0;

{     ********** searchfor{next eigenvalues **********}
60:   if en < low then goto 340;
      its:=0;
      na:=en-1;
      enm2:=na-1;
{     ********** lookfor{single small sub-diagonal element}

70:  for{80} ll:=low to en do
     begin
         
         l:=en+low-ll;
         if l = low then goto 100;
         s:=ABS(h[l-1,l-1])+ABS(h[l,l]);
         if s = 0.0 then s:=norm;
         if ABS(h[l,l-1]) <= (machep*s) then goto 100;
 {80} end;

{     ********** form shift **********}

100:  x:=h[en,en];
       if l = en then goto 270;
      y:=h[na,na];
      w:=h[en,na]*h[na,en];
      if l = na then goto 280;
      if its = 30 then
      begin
        ierr:=en;
        goto 1000;
      end;

{      if (its <> 10) and (its <> 20)  then goto 130;}
      if (its = 10) or (its = 20)  then 
      begin
{     ********** form exceptional shift **********}
        t:=t+x;

        for{120} i:=low to en do
           {120} h[i,i]:=h[i,i]-x;

        s:=ABS(h[en,na])+ABS(h[na,enm2]);
        x:=0.75*s;
        y:=x;
        w:=-0.4375*s*s;
{130:}  its:=its+1;
      end;

{     ********** lookfor{two consecutive small
                 sub-diagonal elements.}

     for{140} mm:=l to enm2 do
     begin
         m:=enm2+l-mm;
         zz:=h[m,m];
         r:=x-zz;
         s:=y-zz;
         p:=(r*s-w)/h[m+1,m]+h[m,m+1];
         q:=h[m+1,m+1]-zz-r-s;
         r:=h[m+2,m+1];
         s:=ABS(p)+ABS(q)+ABS(r);
         p:=p/s;
         q:=q/s;
         r:=r/s;
         if m = l then goto 150;
         if (ABS(h[m,m-1])*(ABS(q)+ABS(r)))
               <= ( machep*ABS(p)
              *(ABS(h[m-1,m-1])+ABS(zz)+ABS(h[m+1,m+1]) )) 
         then goto 150;
 {140} end;

  150: mp2:=m+2;

     for{160} i:=mp2 to en do
     begin
         h[i,i-2]:=0.0;
        { if i = mp2 goto 160;}
         if i <> mp2 then h[i,i-3]:=0.0;
{160}end;

{     ********** double qr step involving rows l to en and
                 columns m to en **********}
     for {260} k:=m to na do
     begin{a}
       notlas:= (k <> na);
        { if k = m goto 170}
       if k <> m then
       begin{b}
         p:=h[k,k-1];
         q:=h[k+1,k-1];
         r:=0.0;
         if notlas then r:=h[k+2,k-1];
         x:=ABS(p)+ABS(q)+ABS(r);
         if x = 0.0 then goto 260;
         p:=p/x;
         q:=q/x;
         r:=r/x;
       end;{b}

       s:=dsign(SQRT(p*p+q*q+r*r),p);
       if k<> m 
          then h[k,k-1]:=-s*x
          else if l <> m then h[k,k-1]:=-h[k,k-1];

       p:=p+s;
       x:=p/s;
       y:=q/s;
       zz:=r/s;
       q:=q/p;
       r:=r/p;

{     ********** row modification **********}

            for{210} j:=k to n do
            begin{d}
              p:=h[k,j]+q*h[k+1,j];
              if notlas then
              begin{e}
                p:=p+r*h[k+2,j];
                h[k+2,j]:=h[k+2,j]-p*zz;
              end;{e}
              h[k+1,j]:=h[k+1,j]-p*y;
              h[k,j]:=h[k,j]-p*x;
     {210}  end;{d}

         j:=min0(en,k+3);

{     ********** column modification **********}

            for{230} i:=1 to j do
            begin{d}
              p:=x*h[i,k]+y*h[i,k+1];
              if notlas then
              begin{e}
                p:=p+zz*h[i,k+2];
                h[i,k+2]:=h[i,k+2]-p*r;
              end;{e}
              h[i,k+1]:=h[i,k+1]-p*q;
              h[i,k]:=h[i,k]-p;
            end;{d}

{     ********** accumulate transformations **********}

          for {250} i:=low to igh do
          begin{d}
            p:=x*z[i,k]+y*z[i,k+1];
            if notlas then
            begin{e}
              p:=p+zz*z[i,k+2];
              z[i,k+2]:=z[i,k+2]-p*r;
            end;{e}
            z[i,k+1]:=z[i,k+1]-p*q;
            z[i,k]:=z[i,k]-p;
   {250}  end;{d}
  260: end;{a}

       goto 70;
{     ********** one root found **********}
 270: h[en,en]:=x+t;
      wr[en]:=h[en,en];
      wi[en]:=0.0;
      en:=na;
   goto 60;
{     ********** two roots found **********}
 280: p:=(y-x)/2.0;
      q:=p*p+w;
      zz:=SQRT(ABS(q));
      h[en,en]:=x+t;
      x:=h[en,en];
      h[na,na]:=y+t;
      if q < 0.0 then goto 320;
{     ********** real pair **********}
      zz:=p+dsign(zz,p);
      wr[na]:=x+zz;
      wr[en]:=wr[na];
      if zz <> 0.0 then wr[en]:=x-w/zz;
      wi[na]:=0.0;
      wi[en]:=0.0;
      x:=h[en,na];
      s:=ABS(x)+ABS(zz);
      p:=x/s;
      q:=zz/s;
      r:=SQRT(p*p+q*q);
      p:=p/r;
      q:=q/r;
{     ********** row modification **********}
     for{290} j:=na to n do
     begin
         zz:=h[na,j];
         h[na,j]:=q*zz+p*h[en,j];
         h[en,j]:=q*h[en,j]-p*zz;
 {290} end;
{     ********** column modification **********}
     for{300} i:=1 to en do
     begin
         zz:=h[i,na];
         h[i,na]:=q*zz+p*h[i,en];
         h[i,en]:=q*h[i,en]-p*zz;
 {300} end;
{     ********** accumulate transformations **********}
     for{310} i:=low to igh do
     begin
         zz:=z[i,na];
         z[i,na]:=q*zz+p*z[i,en];
         z[i,en]:=q*z[i,en]-p*zz;
 {310} end;

   goto 330;
{     ********** complex pair **********}
  320:wr[na]:=x+p;
      wr[en]:=x+p;
      wi[na]:= zz;
      wi[en]:=-zz;
  330:en:=enm2;
   goto 60;
{     ********** all roots found.  backsubstitute to find    }
{                vectors of upper triangular form ********** }
{  340 if norm = 0.0 goto 1001}
  340:if norm <> 0.0 then
    begin{0.5}
     for nn:= 1 to n do
     begin{1}
         en:= n+1-nn;
         p:= wr[en];
         q:= wi[en];
         na:= en-1;
        if q=0 then
        begin{2}{** double vector **}
          m:= en;
          h[en,en]:= 1.0;
          if na<>0 then
          begin{3}
            for ii:= 1 to na do
            begin{4}
              i:= en-ii;

              w:= h[i,i]-p;
              r:= h[i,en];

              if m<=na then
                for j:= m to na do r:= r+h[i,j]*h[j,en];

              if wi[i]<0.0 then
              begin{5}
                zz:= w;
                s:= r;
              end{5} else
              begin{5}
                m:= i;
                if wi[i]=0.0 then
                begin{6}
                  t:= w;
                  if w =0.0 then t:= machep*norm;
                  h[i,en]:= -r/t;
                end{6} else
                begin{6}
     { ********** solve double equations **********}
                  x:= h[i,i+1];
                  y:= h[i+1,i];
                  q:= (wr[i]-p)*(wr[i]-p)+wi[i]*wi[i];
                  t:= (x*s-zz*r)/q;
                  h[i,en]:= t;
                  if Abs(x)<=Abs(zz) 
                     then  h[i+1,en]:= (-s-y*t)/zz
                     else  h[i+1,en]:= (-r-w*t)/x;
                end{6};
              end{5};
            end{4};
        end;{3}
      end{2: of if q=0 then }else
      if q<0 then{ ** Complex ** }
      begin{2}
        m:= na;
{**     ********** last vector component chosen imaginary so that  **}
{**              eigenvector matrix is triangular **********       **}
      { if Abs(h[en,na))<=Abs(h[na,en)) then  720 }
       if na<>0 then
       begin {2.5}        
        if Abs(h[en,na])>Abs(h[na,en]) then 
        begin{3}
          h[na,na]:= q/h[en,na];
          h[na,en]:= -(h[en,en]-p)/h[en,na];
        end{3} else
        begin{3}
          etdiv(z3r,z3i,0.0,-h[na,en],h[na,na]-p,q);
          h[na,na]:= z3r;
          h[na,en]:= z3i;
        end;{3}
        h[en,na]:= 0.0;
        h[en,en]:= 1.0;
        enm2:= na-1;

       {if enm2 =0 then  800}
        if enm2<>0 then 
        begin{3}
         for ii:= 1 to enm2 do
         begin{4}
           i:= na-ii;
           w:= h[i,i]-p;
           ra:= 0.0;
           sa:= h[i,en];

           for j:= m to na do
           begin{5}
             ra:= ra+h[i,j]*h[j,na];
             sa:= sa+h[i,j]*h[j,en];
           end;{5}

           if not (wi[i] >= 0.0) then 
           begin{5}
             zz:= w;
             r:= ra;
             s:= sa;
           end{5} else
           begin{5}
             m:= i;
             if wi[i]=0.0 then
             begin{6}
               etdiv(z3r,z3i,-ra,-sa,w,q);
               h[i,na]:= z3r;
               h[i,en]:= z3i;
             end{6}
             else
{  780     ********** solve complex equations ********** }
             begin{6}
               x:= h[i,i+1];
               y:= h[i+1,i];
               vr:= (wr[i]-p)*(wr[i]-p)+wi[i]*wi[i]-q*q;
               vi:= (wr[i]-p)*2.0*q;

               if (vr =0.0)  and  (vi =0.0)  then
                vr:= machep*norm*(Abs(w)+Abs(q)+Abs(x)+Abs(y)+Abs(zz));

               etdiv(z3r,z3i,x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi);
 
               h[i,na]:= z3r;
               h[i,en]:= z3i;

               if (Abs(x)<=Abs(zz)+Abs(q))  then
               begin{7}
                 etdiv(z3r,z3i,-r-y*h[i,na],-s-y*h[i,en],zz,q);
                 h[i+1,na]:= z3r;
                 h[i+1,en]:= z3i;
               end{7} else
               begin{7}
                 h[i+1,na]:= (-ra-w*h[i,na]+q*h[i,en])/x;
                 h[i+1,en]:= (-sa-w*h[i,en]-q*h[i,na])/x;
               end;{7}
             end;{6}
           end;{4}
         end;{5}
       end{3};
      end;{2.5}
     end;{2}
   end;{1}


{*     ********** end back substitution.               *}
{*               vectors of isolated roots **********  *}
      for i:= 1 to n do
         if not((i>=low) and (i<=igh))  then  
            for j:= i to n do z[i,j]:= h[i,j];

{*     ********** multiply by transformation matrix to give  *}
{*                vectors of original full matrix.           *}
{*                for j:=n step -1 until low for --          *}


      for jj:= low to n do
      begin{1}
         j:= n+low-jj;
         m:= min0(j,igh);
         for i:= low to igh do
         begin{2}
            zz:= 0.0;
            for k:= low to  m do zz:= zz+z[i,k]*h[k,j];
            z[i,j]:= zz;
         end;{2}
      end;{1}
    end;{0.5}
    en:=0;
{     ********** last card of ety2 **********}
  1000:ierr:=en;
 end;

  end.
