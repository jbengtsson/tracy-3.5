	subroutine graheds(str, len, h)

	implicit none

	character	str*80, chatra*80
	real		h
	integer		len, i

        chatra=' '
	do 20 i=1,len
          chatra(i:i)=str(i:i)
20	continue
	call grahed(chatra, h)
        return
        end


	subroutine gracoms(h, xstr, xlen, ixpos, ystr, ylen, nyelm)

	implicit none

	character	xstr*80, ystr*80, fxstr*80, fystr*80
	real		h
	integer		xlen, ylen, ixpos, nyelm, i

        fxstr=' '
	fystr=' '
	do 20 i=1,xlen
          fxstr(i:i)=xstr(i:i)
20	continue
	do 40 i=1,ylen
          fystr(i:i)=ystr(i:i)
40	continue
	call gracom(h, fxstr, ixpos, fystr, 1)
        return
        end


	subroutine gratexs(xst, yst, h, str, len, win, xe, ye, jfplt)

	implicit none

	character	str*80, chatra*80
	real		xst, yst, h, win, xe, ye
	integer		len, jfplt, i

        chatra=' '
	do 20 i=1,len
          chatra(i:i)=str(i:i)
20	continue
	call gratex(xst, yst, h, chatra, win, xe, ye, jfplt)
        return
        end
