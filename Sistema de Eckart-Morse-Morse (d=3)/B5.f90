program CNF3DoF

implicit none

real(8), allocatable		::	Hqp(:,:),temp_alloc(:,:),WnQP(:,:)
complex(8), allocatable		::	Hxxi(:,:),Hn(:,:),Wn(:,:)
real(8)				::	moldeQP(12174,8),temp_real(12174,8),temp_moyal(12174,8),pol_id(1,8),t1
complex(8)			::	moldeXXI(12174,8),temp_cplx(12174,8)
real(8), parameter		::	lambda=0.7349552361081487d0,den3=-1.d0/(24.d0),den5=1.d0/(16.d0*120.d0)
complex(8), parameter		::	iw2=dcmplx(0.d0,1.267290444967991d0),iw3=dcmplx(0.d0,1.822517936036739d0)
integer				::	i1,i2,i3,i4,i5,i6,i7,k,ind(0:5,0:10,0:10,0:10,0:10,0:10,0:10)=0,n,cont,s,l,c3=0,c5=0, &
					ex3(108,7),ex5(1458,7),ord
integer, parameter		::	m=10

print *
print *, "Inicializando..."

! leer polinomio inicial

allocate (Hqp(8008,8))
 open (unit=99, file='data.txt', status='old', action='read')
do k=1,8008
	read(99,*) Hqp(k,:)
end do
 close(unit=99, status='KEEP') 

! definir polinomios "molde" en coordenadas (q,p) & (x,xi)

k=0
do i1=0,5
	do i2=0,10-2*i1
		do i3=0,10-2*i1-i2
			do i4=0,10-2*i1-i2-i3
				do i5=0,10-2*i1-i2-i3-i4
					do i6=0,10-2*i1-i2-i3-i4-i5
						do i7=0,10-2*i1-i2-i3-i4-i5-i6
							k = k + 1
							ind(i1,i2,i3,i4,i5,i6,i7)=k
							moldeQP(k,:) = (/ dble(i1),dble(i2),dble(i3),dble(i4),dble(i5),dble(i6),dble(i7),0.d0 /)
							moldeXXI(k,:) = (/ dcmplx(i1),dcmplx(i2),dcmplx(i3),dcmplx(i4),dcmplx(i5),dcmplx(i6),dcmplx(i7),dcmplx(0) /)
						end do
					end do
				end do
			end do
		end do
	end do
end do

! hallar coef de derivadas a utilizar en el cálculos de los corchetes de Moyal

do i1=1,3
	do i2=1,3
		do i3=1,3
			do k=0,3
				c3=c3+1
				ex3(c3,1:6) = (/0,0,0,0,0,0/)
				ex3(c3,7) = k
				if (k==0) then
				ex3(c3,i1) = ex3(c3,i1) + 1
				ex3(c3,i2) = ex3(c3,i2) + 1
				ex3(c3,i3) = ex3(c3,i3) + 1
				else if (k==1) then
				ex3(c3,i1) = ex3(c3,i1) + 1
				ex3(c3,i2) = ex3(c3,i2) + 1
				ex3(c3,i3+3) = ex3(c3,i3+3) + 1
				else if (k==2) then
				ex3(c3,i1) = ex3(c3,i1) + 1
				ex3(c3,i2+3) = ex3(c3,i2+3) + 1
				ex3(c3,i3+3) = ex3(c3,i3+3) + 1
				else
				ex3(c3,i1+3) = ex3(c3,i1+3) + 1
				ex3(c3,i2+3) = ex3(c3,i2+3) + 1
				ex3(c3,i3+3) = ex3(c3,i3+3) + 1
				endif
			enddo
		enddo
	enddo
enddo

do i1=1,3
	do i2=1,3
		do i3=1,3
			do i4=1,3
				do i5=1,3
					do k=0,5
						c5=c5+1
						ex5(c5,1:6) = (/0,0,0,0,0,0/)
						ex5(c5,7) = k
						if (k==0) then
							ex5(c5,i1) = ex5(c5,i1) + 1
							ex5(c5,i2) = ex5(c5,i2) + 1
							ex5(c5,i3) = ex5(c5,i3) + 1
							ex5(c5,i4) = ex5(c5,i4) + 1
							ex5(c5,i5) = ex5(c5,i5) + 1
						else if (k==1) then
							ex5(c5,i1) = ex5(c5,i1) + 1
							ex5(c5,i2) = ex5(c5,i2) + 1
							ex5(c5,i3) = ex5(c5,i3) + 1
							ex5(c5,i4) = ex5(c5,i4) + 1
							ex5(c5,i5+3) = ex5(c5,i5+3) + 1
						else if (k==2) then
							ex5(c5,i1) = ex5(c5,i1) + 1
							ex5(c5,i2) = ex5(c5,i2) + 1
							ex5(c5,i3) = ex5(c5,i3) + 1
							ex5(c5,i4+3) = ex5(c5,i4+3) + 1
							ex5(c5,i5+3) = ex5(c5,i5+3) + 1
						else if (k==3) then
							ex5(c5,i1) = ex5(c5,i1) + 1
							ex5(c5,i2) = ex5(c5,i2) + 1
							ex5(c5,i3+3) = ex5(c5,i3+3) + 1
							ex5(c5,i4+3) = ex5(c5,i4+3) + 1
							ex5(c5,i5+3) = ex5(c5,i5+3) + 1
						else if (k==4) then
							ex5(c5,i1) = ex5(c5,i1) + 1
							ex5(c5,i2+3) = ex5(c5,i2+3) + 1
							ex5(c5,i3+3) = ex5(c5,i3+3) + 1
							ex5(c5,i4+3) = ex5(c5,i4+3) + 1
							ex5(c5,i5+3) = ex5(c5,i5+3) + 1
						else
							ex5(c5,i1+3) = ex5(c5,i1+3) + 1
							ex5(c5,i2+3) = ex5(c5,i2+3) + 1
							ex5(c5,i3+3) = ex5(c5,i3+3) + 1
							ex5(c5,i4+3) = ex5(c5,i4+3) + 1
							ex5(c5,i5+3) = ex5(c5,i5+3) + 1
						endif
					enddo
				enddo
			enddo
		enddo
	enddo
enddo

! definir el polinomio identidad (p(x,xi) = 1)

pol_id(1,7) = moldeQP(1,7)
pol_id(1,8) = 1.d0

! iniciar procedimiento iterativo

DO n=3,m


	print *
	print *, '<--------------> ',n,'           <-------------->'
	print *
	print *, 'Cambiando coordenadas: Hqp -> Hxxi'

	call QPtoXXI(Hqp,temp_cplx)

	if (allocated(Hxxi) .eqv. .true.) then
		deallocate (Hxxi)
	end if

	call reduceXXI(temp_cplx,Hxxi)

	print *, 'Obteniendo Hn'

	if (allocated(Hn) .eqv. .true.) then
		deallocate (Hn)
	end if

	call OrdenXXI(Hxxi,n,Hn)

	print *, 'Obteniendo Wn'

	cont = 0
	do k=1,size(Hn(:,1))
		if (   nint(real(Hn(k,2))).ne.nint(real(Hn(k,5))) .OR. &
		       nint(real(Hn(k,3))).ne.nint(real(Hn(k,6))) .OR. &
		       nint(real(Hn(k,4))).ne.nint(real(Hn(k,7)))   ) then
			cont = cont + 1
			temp_cplx(cont,1:7) = Hn(k,1:7)
			temp_cplx(cont,8) = Hn(k,8)/( lambda*(nint(real(Hn(k,5)))-nint(real(Hn(k,2)))) &
			+ iw2*(nint(real(Hn(k,6)))-nint(real(Hn(k,3)))) &
			+ iw3*(nint(real(Hn(k,7)))-nint(real(Hn(k,4)))) )
		end if
	end do

	if (allocated(Wn) .eqv. .true.) then
		deallocate (Wn)
	end if

	call reduceXXI(temp_cplx(1:cont,:),Wn)

! expresar Wn en coordenadas (q,p)

	call XXItoQP(Wn,temp_real)

	if (allocated(WnQP) .eqv. .true.) then
		deallocate (WnQP)
	end if

	call reduceQP(temp_real,WnQP)

! determinar forma normal de orden n

	print *, 'Calculando forma normal de orden',n

	temp_real(:,:) = moldeQP(:,:)

	do s=n,m

		print *, 's =',s

		if (allocated(temp_alloc) .eqv. .true.) then
			deallocate (temp_alloc)
		end if

		call OrdenQP(Hqp,s,temp_alloc)
		call pol_prod(temp_alloc,pol_id,temp_real)

		do k=1,floor(real(s)/real(n-2))
			print *, 'k =',k
			deallocate (temp_alloc)
			ord = s-k*(n-2)
			call OrdenQP(Hqp,ord,temp_alloc)

			do l=1,k
				call MoyalBracket(WnQP,temp_alloc,MIN(n,ord),temp_moyal)
				ord = ord + n-2
				deallocate (temp_alloc)
				call reduceQP(temp_moyal,temp_alloc)
			end do

		call ct_prod(temp_moyal,1/dble(factorial(k)))
		call pol_sum(temp_real,temp_moyal)

		end do

	end do

! modificar los términos de Hqp de orden >= n

	do k=1,size(Hqp(:,1))
		if (nint(2*Hqp(k,1)+Hqp(k,2)+Hqp(k,3)+Hqp(k,4)+Hqp(k,5)+Hqp(k,6)+Hqp(k,7)) >= n) then
			Hqp(k,8) = 0.d0
		end if
	end do

	deallocate (temp_alloc)
	call reduceQP(Hqp,temp_alloc)
	temp_moyal(:,:) = moldeQP(:,:)
	call pol_prod(temp_alloc,pol_id,temp_moyal)

	call pol_sum(temp_real,temp_moyal)

	deallocate (Hqp)
	call reduceQP(temp_real,Hqp)

	print *, 'Fin de iteración'
END DO



call QPtoXXI(Hqp,temp_cplx)
if (allocated(Hxxi) .eqv. .true.) then
	deallocate (Hxxi)
end if
call reduceXXI(temp_cplx,Hxxi)

do k=0,m

	print *
	print *, '--------------------------------------------------------------------------------------------'
	print *
	
	if (allocated(Hn) .eqv. .true.) then
		deallocate (Hn)
	end if

	call OrdenXXI(Hxxi,k,Hn)

	do l=1,size(Hn(:,1))
		print *, Hn(l,:)
	end do

end do

 open (unit=11, file="--output--QNF.txt")
do k=1,size(Hxxi(:,1))
	write(11, '(*(F20.16 : ", "))') Hxxi(k,:)
enddo
 close(unit=11, status='KEEP') 

contains
										! CALCULAR DERIVADAS
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

subroutine derivar(pol,e1,e2,e3,e4,e5,e6,out_pol)

real(8), allocatable, intent(in)	::	pol(:,:)
integer, intent(in)			::	e1,e2,e3,e4,e5,e6
real(8)					::	temp
integer 				::	k,i
real(8), allocatable, intent(out)	::	out_pol(:,:)

allocate (out_pol(size(pol(:,1)),8))

do k=1,size(pol(:,1))

	if (nint(pol(k,2)) >= e1 .AND. nint(pol(k,3)) >= e2 .AND. nint(pol(k,4)) >= e3 .AND. &
	   nint(pol(k,5)) >= e4 .AND. nint(pol(k,6)) >= e5 .AND. nint(pol(k,7)) >= e6) then
   
		i=ind(nint(pol(k,1)),nint(pol(k,2)) - e1,nint(pol(k,3)) - e2,&
		   nint(pol(k,4)) - e3,nint(pol(k,5)) - e4,nint(pol(k,6)) - e5,nint(pol(k,7)) - e6)
		out_pol(k,1:7) = moldeQP(i,1:7)
		temp = 1.d0

		do i=nint(pol(k,2))+1-e1,nint(pol(k,2))
			temp = temp*i
		end do
		do i=nint(pol(k,3))+1-e2,nint(pol(k,3))
			temp = temp*i
		end do
		do i=nint(pol(k,4))+1-e3,nint(pol(k,4))
			temp = temp*i
		end do
		do i=nint(pol(k,5))+1-e4,nint(pol(k,5))
			temp = temp*i
		end do
		do i=nint(pol(k,6))+1-e5,nint(pol(k,6))
			temp = temp*i
		end do
		do i=nint(pol(k,7))+1-e6,nint(pol(k,7))
			temp = temp*i
		end do

		out_pol(k,8) = pol(k,8)*temp
	else
		out_pol(k,1:7) = pol(k,1:7)
		out_pol(k,8) = 0.d0
	end if

end do

end subroutine
										! MULTIPLICAR POLINOMIOS
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

subroutine pol_prod(polA,polB,prev_sum)

real(8), intent(in)		::	polA(:,:),polB(:,:)
real(8), intent(inout)		::	prev_sum(12174,8)
integer				::	i,j,k

do i=1,size(polA(:,1))
	do j=1,size(polB(:,1))
		if (nint(2.d0*(polA(i,1)+polB(j,1))+(polA(i,2)+polB(j,2))+(polA(i,3)+polB(j,3))+(polA(i,4)+polB(j,4))+ &
		  	(polA(i,5)+polB(j,5))+(polA(i,6)+polB(j,6))+(polA(i,7)+polB(j,7))) .LE. 10 ) then
			k=ind(nint(polA(i,1)+polB(j,1)),nint(polA(i,2)+polB(j,2)),nint(polA(i,3)+polB(j,3)),nint(polA(i,4)+polB(j,4)), &
		  	nint(polA(i,5)+polB(j,5)),nint(polA(i,6)+polB(j,6)),nint(polA(i,7)+polB(j,7)))
			prev_sum(k,8) = prev_sum(k,8) + polA(i,8)*polB(j,8)
		end if
	end do
end do

end subroutine
										! SUMAR POLINOMIOS
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

subroutine pol_sum(polA,polB)

real(8), intent(in)		::	polB(12174,8)
real(8), intent(inout)		::	polA(12174,8)
integer				::	k

polA(:,8) = polA(:,8) + polB(:,8)

end subroutine
										! MULTIPLICAR POR CONSTANTES
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

subroutine ct_prod(pol,cte)

real(8), intent(inout)	::	pol(:,:)
real(8), intent(in)	::	cte

pol(:,8) = pol(:,8)*cte

end subroutine
										! ELEGIR MONOMIOS DE ORDEN N, coordenadas (q,p)
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

subroutine OrdenQP(pol,ord,out_pol)

real(8), allocatable, intent(in)	::	pol(:,:)
real(8), allocatable, intent(out)	::	out_pol(:,:)
integer, intent(in)			::	ord
integer					::	k,cont
cont=0
do k=1,size(pol(:,1)) 
	if (nint(real(2*pol(k,1)+pol(k,2)+pol(k,3)+pol(k,4)+pol(k,5)+pol(k,6)+pol(k,7))) == ord) then
		cont = cont + 1
	end if
end do

if (cont > 0) then

	allocate (out_pol(cont,8))

	cont = 0
	do k=1,size(pol(:,1)) 
		if (nint(real(2*pol(k,1)+pol(k,2)+pol(k,3)+pol(k,4)+pol(k,5)+pol(k,6)+pol(k,7))) == ord) then
			cont = cont + 1
			out_pol(cont,:) = pol(k,:)
		end if
	end do

else
	allocate (out_pol(1,8))
	out_pol(1,1) = 0.d0
	out_pol(1,2) = 0.d0
	out_pol(1,3) = 0.d0
	out_pol(1,4) = 0.d0
	out_pol(1,5) = 0.d0
	out_pol(1,6) = 0.d0
	out_pol(1,7) = 0.d0
	out_pol(1,8) = 0.d0
end if

end subroutine
										! ELEGIR MONOMIOS DE ORDEN N, coordenadas (x,xi)
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

subroutine OrdenXXI(pol,ord,out_pol)

complex(8), allocatable, intent(in)	::	pol(:,:)
complex(8), allocatable, intent(out)	::	out_pol(:,:)
integer, intent(in)			::	ord
integer					::	k,cont
cont=0
do k=1,size(pol(:,1)) 
	if (nint(real(2*pol(k,1)+pol(k,2)+pol(k,3)+pol(k,4)+pol(k,5)+pol(k,6)+pol(k,7))) == ord) then
		cont = cont + 1
	end if
end do

if (cont > 0) then

	allocate (out_pol(cont,8))

	cont = 0
	do k=1,size(pol(:,1)) 
		if (nint(real(2*pol(k,1)+pol(k,2)+pol(k,3)+pol(k,4)+pol(k,5)+pol(k,6)+pol(k,7))) == ord) then
			cont = cont + 1
			out_pol(cont,:) = pol(k,:)
		end if
	end do

else
	allocate (out_pol(1,8))
	out_pol(1,1) = dcmplx(0,0)
	out_pol(1,2) = dcmplx(0,0)
	out_pol(1,3) = dcmplx(0,0)
	out_pol(1,4) = dcmplx(0,0)
	out_pol(1,5) = dcmplx(0,0)
	out_pol(1,6) = dcmplx(0,0)
	out_pol(1,7) = dcmplx(0,0)
	out_pol(1,8) = dcmplx(0,0)
end if

end subroutine
										! SIMPLIFICAR POLINOMIOS, coordenadas (q,p)
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

subroutine reduceQP(pol,out_pol)

real(8), intent(in)				::	pol(:,:)
real(8), parameter				::	tol=1.d-12
real(8), allocatable, intent(out)		::	out_pol(:,:)
integer						::	k,cont
cont=0
do k=1,size(pol(:,1))
	if (abs(pol(k,8)) > tol) then
		cont = cont + 1
	end if
end do

if (cont > 0) then
	allocate (out_pol(cont,8))
	cont = 0
	do k=1,size(pol(:,1))
		if (abs(pol(k,8)) > tol) then
			cont = cont + 1
			out_pol(cont,:) = pol(k,:)
		end if
	end do
else
	allocate (out_pol(1,8))
	out_pol(1,1) = 0.d0
	out_pol(1,2) = 0.d0
	out_pol(1,3) = 0.d0
	out_pol(1,4) = 0.d0
	out_pol(1,5) = 0.d0
	out_pol(1,6) = 0.d0
	out_pol(1,7) = 0.d0
	out_pol(1,8) = 0.d0
end if

end subroutine
										! SIMPLIFICAR POLINOMIOS, coordenadas (x,xi)
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

subroutine reduceXXI(pol,out_pol)

complex(8), intent(in)				::	pol(:,:)
real(8), parameter				::	tol=1.d-12
complex(8), allocatable, intent(out)		::	out_pol(:,:)
integer						::	k,cont
cont = 0
do k=1,size(pol(:,1))
	if (abs(pol(k,8)) > tol) then
		cont = cont + 1
	end if
end do

if (cont > 0) then
	allocate (out_pol(cont,8))
	cont = 0
	do k=1,size(pol(:,1))
		if (abs(pol(k,8)) > tol) then
			cont = cont + 1
			out_pol(cont,:) = pol(k,:)
		end if
	end do
else
	allocate (out_pol(1,8))
	out_pol(1,1) = 0.d0
	out_pol(1,2) = 0.d0
	out_pol(1,3) = 0.d0
	out_pol(1,4) = 0.d0
	out_pol(1,5) = 0.d0
	out_pol(1,6) = 0.d0
	out_pol(1,7) = 0.d0
	out_pol(1,8) = 0.d0
end if

end subroutine
										! DEFINICIONES: FACTORIAL, COEF BINOMIALES
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

  function factorial(n) result (res)
 
    implicit none
    integer, intent (in) :: n
    integer :: res
    integer :: i
 
    res = product ((/(i, i = 1, n)/))
 
  end function factorial
 
  function choose(n, k) result (res)
 
    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: k
    integer :: res
 
    res = factorial (n) / (factorial (k) * factorial (n - k))
 
  end function choose
										! CAMBIO DE COORDENADAS: (q,p) -> (x,xi)
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

subroutine QPtoXXI(pol,out_pol)

real(8), allocatable, intent(in)	::	pol(:,:)
complex(8), intent(out)			::	out_pol(12174,8)
integer					::	k,m,e2,e3,e5,e6

out_pol(:,:) = moldeXXI(:,:)

do k=1,size(pol(:,1))
	do e2=0,nint(pol(k,3))
		do e3=0,nint(pol(k,4))
			do e5=0,nint(pol(k,6))
				do e6=0,nint(pol(k,7))
					m=ind(nint(pol(k,1)),nint(pol(k,2)),e2+nint(pol(k,6))-e5,&
					  e3+nint(pol(k,7))-e6,nint(pol(k,5)),nint(pol(k,3))-e2+e5,nint(pol(k,4))-e3+e6)
					out_pol(m,8) = out_pol(m,8)+pol(k,8)*(1.d0/2.d0)**(nint((pol(k,3))+nint(pol(k,4))+ &
					nint(pol(k,6))+nint(pol(k,7)))/2.d0)* &
					choose(nint(pol(k,3)),e2)*choose(nint(pol(k,4)),e3)*choose(nint(pol(k,6)),e5)* &
					choose(nint(pol(k,7)),e6)* &
					dcmplx(0,1)**(nint(pol(k,3))+nint(pol(k,4))-e2-e3+nint(pol(k,6))-e5+nint(pol(k,7))-e6)
				end do
			end do
		end do
	end do
end do

end subroutine
										! CAMBIO DE COORDENADAS: (x,xi) -> (q,p)
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

subroutine XXItoQP(pol,out_pol)

complex(8), intent(in)		::	pol(:,:)
real(8), intent(out)		::	out_pol(12174,8)
integer				::	k,m,e2,e3,e5,e6

out_pol(:,:) = moldeQP(:,:)

do k=1,size(pol(:,1))
	do e2=0,nint(real(pol(k,3)))
		do e3=0,nint(real(pol(k,4)))
			do e5=0,nint(real(pol(k,6)))
				do e6=0,nint(real(pol(k,7)))
					m=ind(nint(real(pol(k,1))),nint(real(pol(k,2))),e2+nint(real(pol(k,6)))-e5,&
					  e3+nint(real(pol(k,7)))-e6,nint(real(pol(k,5))),nint(real(pol(k,3)))-e2+e5,nint(real(pol(k,4)))-e3+e6)
					out_pol(m,8) = out_pol(m,8)+pol(k,8)*(1.d0/2.d0)**(nint(real((pol(k,3)))+ &
					nint(real(pol(k,4)))+ &
					nint(real(pol(k,6)))+nint(real(pol(k,7))))/2.d0)* &
					choose(nint(real(pol(k,3))),e2)*choose(nint(real(pol(k,4))),e3)* &
					choose(nint(real(pol(k,6))),e5)* &
					choose(nint(real(pol(k,7))),e6)* &
					dcmplx(0,-1)**(nint(real(pol(k,3)))+nint(real(pol(k,4)))- &
					e2-e3+nint(real(pol(k,6)))-e5+nint(real(pol(k,7)))-e6)					
				end do
			end do
		end do
	end do
end do

end subroutine
										! CORCHETE DE MOYAL
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

subroutine MoyalBracket(polA,polB,ord,out_pol)

real(8), allocatable, intent(in)	::	polA(:,:),polB(:,:)
real(8), intent(out)			::	out_pol(12174,8)
real(8), allocatable			::	dz1(:,:),dz2(:,:),temp(:,:),mono(:,:)
integer, parameter			::	zeros(6)=(/0,0,0,0,0,0/)
integer, intent(in)			::	ord
integer					::	k,l,ex(6)

allocate (temp(12174,8))
temp(:,:) = moldeQP(:,:)
out_pol(:,:) = moldeQP(:,:)
ex(:) = zeros(:)

! corchete de Poisson

do k=1,6

	if (allocated(dz1) .eqv. .true.) then
		deallocate (dz1)
		deallocate (dz2)
	end if

	ex(k) = 1
	call derivar(polA,ex(1),ex(2),ex(3),ex(4),ex(5),ex(6),dz1)
	ex(k) = 0

	if (k<4) then
		ex(k+3) = 1
		call derivar(polB,ex(1),ex(2),ex(3),ex(4),ex(5),ex(6),dz2)
		ex(k+3) = 0
	else
		ex(k-3) = 1
		call derivar(polB,ex(1),ex(2),ex(3),ex(4),ex(5),ex(6),dz2)
		call ct_prod(dz1,-1.d0)
		ex(k-3) = 0

	end if

	call pol_prod(dz1,dz2,out_pol)	

end do

! derivadas de tercer orden

if (ord >= 3) then

	do k=1,108
	
		deallocate (dz1)
		deallocate (dz2)

		ex(:) = ex3(k,1:6)
		call derivar(polA,ex(1),ex(2),ex(3),ex(4),ex(5),ex(6),dz1)
		call ct_prod(dz1,(-1.d0)**ex3(k,7)*dble(choose(3,ex3(k,7))))
		ex(1:3) = ex3(k,4:6)
		ex(4:6) = ex3(k,1:3)
		call derivar(polB,ex(1),ex(2),ex(3),ex(4),ex(5),ex(6),dz2)

		call pol_prod(dz1,dz2,temp)

	end do

	allocate (mono(1,8))
	mono(1,:) = (/ 2.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,den3 /)

	call pol_prod(mono,temp,out_pol)

end if

! derivadas de quinto orden

if (ord >= 5) then

	temp(:,:) = moldeQP(:,:)
	do k=1,1458
	
		deallocate (dz1)
		deallocate (dz2)

		ex(:) = ex5(k,1:6)
		call derivar(polA,ex(1),ex(2),ex(3),ex(4),ex(5),ex(6),dz1)
		call ct_prod(dz1,(-1.d0)**ex5(k,7)*dble(choose(5,ex5(k,7))))
		ex(1:3) = ex5(k,4:6)
		ex(4:6) = ex5(k,1:3)
		call derivar(polB,ex(1),ex(2),ex(3),ex(4),ex(5),ex(6),dz2)

		call pol_prod(dz1,dz2,temp)

	end do

	deallocate (mono)
	allocate (mono(1,8))
	mono(1,:) = (/ 4.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,den5 /)

	call pol_prod(mono,temp,out_pol)

end if

end subroutine

! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !
! #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# !

end program
