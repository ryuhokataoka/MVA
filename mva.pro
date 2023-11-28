PRO mva, b1,b2,b3, n1,n2,n3, l1,l2,l3

;+
; NAME: MVA
;
; PURPOSE:
;	This procedure returns the normal vector
;	derived by minimum variance analysis (MVA)
;	of Sonnerup and Carhill, JGR, 171-183, 1967.
;
; INPUT:
;	b1,2,3:	Magnetic field (bx,by,bz) in Cartesian coordinates.
;
; OUTPUT:
;	n1,2,3:	Determined normal vectors (max,med,min).
;	l1,2,3: Determined eigenvalues (max,med,min).
;
; USAGE:
;	MVA_normal = n3
;	MVA_ratio = l2/l3
;
;	by Ryuho Kataoka, Feb 2, 2003
;-

nx=n_elements(b1)
ny=n_elements(b2)
nz=n_elements(b3)

if nx ne ny then stop	;"should have the same dimentions"
if nx ne nz then stop
if ny ne nz then stop

n=nx

m11=total(b1*b1)/n - total(b1)/n*total(b1)/n
m12=total(b1*b2)/n - total(b1)/n*total(b2)/n
m13=total(b1*b3)/n - total(b1)/n*total(b3)/n
m22=total(b2*b2)/n - total(b2)/n*total(b2)/n
m23=total(b2*b3)/n - total(b2)/n*total(b3)/n
m33=total(b3*b3)/n - total(b3)/n*total(b3)/n

m21=m12
m31=m13
m32=m23

M=[	[ m11, m12, m13],$
	[ m21, m22, m23],$
	[ m31, m32, m33]]

A=M

trired, A, D, E	;compute the tridiagonal form of A
triql, D, E, A	;compute the eigenvalues(D) and the eigenvectors(A)

indx=reverse(sort(D))

n1=A(*,indx(0));max
n2=A(*,indx(1));med
n3=A(*,indx(2));min

l1=D(indx(0));max
l2=D(indx(1));med
l3=D(indx(2));min

return
END
