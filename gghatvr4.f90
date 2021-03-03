!ifort -r8  gghatvr3.f90  -mkl -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -o gghatvr3
!ifort -r8 -qopenmp gghatvr3.f90  -mkl -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -o gghatvr3

!include 'clusterm.f90'
program main
!use clusterm
USE lapack95, ONLY: potrf, potri
use blas95, only:  dot

integer, allocatable :: nfrq(:),imean(:),id(:),iifix(:,:,:),iorder(:)
real, allocatable :: y(:),wt(:),wt2(:),sol(:),afix(:,:),frq(:),Gmat(:,:),ifix(:,:,:)
real*4, allocatable :: Amat(:,:)
real :: add_diag(1)
integer :: pnloc(100),pnmeans,pnhap,pndat
character(len=100) :: genotypefile,file_name,option,fil_format,mapfile
character(len=5000000) :: line,char*1,line2*200
character(len=1), allocatable :: allele1(:)
integer :: ids(1000000),idinv(100000),nloc_chr(100)
integer*1, allocatable :: ifour(:)
integer*4 :: jpos(4)=(/0,2,4,6/),jval(4)
logical :: calc_freqs !plinkbed_file, ldmipmrk_file, ldmipout_fil
namelist/inputs/nchroms,pnloc,pndat,genotypefile,mapfile,fil_format,igiv,add_diag,option,ivanraden,calc_freqs

!read steer file
print *,' reading gghatvr4.inp'
open(2,file='gghatvr4.inp',status='old')
nchroms=0; pnloc=0; pndat=0; genotypefile="0"; mapfile="0";fil_format="";add_diag=0.0;option="";ivanraden=1;calc_freqs=.true.
read(2,nml=inputs)
igiv=0        
if(option(1:9)=="noGnoGINV")then
   igiv=-1
elseif(option(1:5)=="Gonly")then
   igiv=0         
elseif(option(1:8)=="GandGINV")then
   igiv=1
endif
if(fil_format(1:5)=="plink")then
   nbytes=ceiling(pdat/4.)
   allocate(ifour(nbytes))
endif   

!pnloc=0
!do i=1,nchroms
!  read(2,*)pnloc(i)
!enddo
!read(2,*)pndat
!
!read(2,'(a)')filnam  !ldmip3.out files. filnam contains a %-sign which will be replaced by the chrom number
!fil_format=""
!IF(FILNAM(1:1)/="0")THEN
!ic=scan(filnam," ")
!if(ic>1)filnam=filnam(1:ic)
!if(index(filnam,'ldmip')>0)then
!  if(index(filnam,'.out')>0)option="ldmipout"; print *,'ldmip.out inputfile'
!endif
!if(index(filnam,'.mrk')>0)then
!  option="ldmipmrk"; print *,'ldmip.mrk inputfile'
!elseif(index(filnam,'.gen')>0)then 
!  option="ldmipgen"; ioption=index(filnam,'.gen'); print *,'ldmip.gen inputfile'
!elseif(index(filnam,'.012')>0)then
!  option="012"; ioption=index(filnam,'.012'); print *,'*.012 inputfile'
!elseif(index(filnam,'.phased')>0)then
!  option="bgl"; ioption=index(filnam,'.bgl'); print *,'beagle format .phased file'
!endif
!ENDIF
!read(2,*)filnam2
!ic=scan(filnam2," ")
!if(ic>1)filnam2=filnam2(1:ic)
ifilnam2=0;
if(mapfile(1:1)/="0" .and. fil_format=="ldmipgen")ifilnam2=1
!read(2,*)igiv,xxxx  !igiv=0/1 indicates whether inverse is needed or not
!if(igiv==1)add_diag(1)=xxxx
!!if(igiv==-2)r2crit=xxxx
!read(2,*)ivanraden   !VanRaden method : 1 or 2 [or -1 or -2 where minus indicates that all allele frequencies are read from file gghatvr3.freq]
!if(igiv==-1)ivanraden=sign(2,ivanraden)
if(.not.(calc_freqs))open(80,file="gghatvr3.freq",status="old")
IF(NCHROMS==0)GOTO 23

ic=scan(genotypefile,' ')
if(ic>1)genotypefile(ic:)=''  !remove perhaps some text at the end
ic=scan(genotypefile,'%')
close(2)

  allocate(y(pndat),wt(pndat*2),afix(maxval(pnloc),pndat),frq(maxval(pnloc)),Gmat(pndat*2,pndat*2))
  allocate(Amat(pndat*2,pndat*2),wt2(pndat*2),allele1(maxval(pnloc)),iorder(maxval(pnloc)))
  Gmat=0.0; Amat=0.; 
  allocate(ifix(2,maxval(pnloc),pndat),iifix(2,maxval(pnloc),pndat),nfrq(maxval(pnloc)),imean(pndat),id(pndat*2))
!set up Amat

IF(fil_format(1:8)=="ldmipout" .and. igiv>=0)then
  open(19,file='ldmip3.ped',status='old')
  do i=1,pndat
    read(19,*)ii,iis,iid
    iis=max(iis,0); iid=max(iid,0)
    i2=i*2; i1=i2-1; iis2=iis*2; iis1=iis2-1; iid2=iid*2; iid1=iid2-1
    Gmat(i1,i1)=1.; Gmat(i2,i2)=1.
    if(iis>0)then
      Gmat(i1,1:i1-1)=.5*Gmat(iis1,1:i1-1)+.5*Gmat(iis2,1:i1-1)
      Gmat(1:i1-1,i1)=Gmat(i1,1:i1-1)
      Gmat(i2,1:i2-1)=.5*Gmat(iid1,1:i2-1)+.5*Gmat(iid2,1:i2-1)
      Gmat(1:i2-1,i2)=Gmat(i2,1:i2-1)
    endif
  enddo
  close(19)
  id(1:pndat)=(/(i,i=1,pndat*2,2)/)
  id(pndat+1:2*pndat)=(/(i,i=2,pndat*2,2)/)
  forall(ii=1:pndat*2,jj=1:pndat*2)Amat(ii,jj)=Gmat(id(ii),id(jj)) !Amat(1:2*pndat,1:2*pndat)=Amat(id(1:2*pndat),id(1:2*pndat))
  Gmat=0.0  
ENDIF !fil_format

Gmat=0
denomsum=0.0
nnloc=0; 
DO ICHROM=1,nchroms
nloc=pnloc(ichrom); ilocus=0
print *,' nloc=',nloc,' no of recs',pndat


if(ic==0)then
  file_name=trim(genotypefile)
elseif(ichrom<=9)then
  file_name=genotypefile
  write(file_name(ic:ic),'(i1)')ichrom
else
  file_name=genotypefile
  file_name(ic+1:)=genotypefile(ic:)
  write(file_name(ic:ic+1),'(i2)')ichrom
endif
if(fil_format(1:5)=="plink")then
   open(4,file=genotypefile,access="stream",form='unformatted',status='old') !read plink.bed file
   read(4)ifour(1:3) !read 3 magic plink numbers
   file_name=genotypefile
   icc=scan(file_name,'.',back=.true.)
   if(icc>0)then
      file_name(icc+1:icc+3)="fam"
      open(5,file=file_name,status="old",iostat=ierr)  !read plink.fam file
      if(ierr==0)then
         print *,"reading ",file_name 
         do i=1,pndat
            read(5,*)idum,ids(i)
         enddo
      else
         ids(1:pndat)=(/(i,i=1,pndat)/)
      endif
   endif
else   
   open(4,file=file_name,status='old') !ldmip3.out
endif   
if(fil_format(1:3)=="bgl")read(4,*)  !skip header line

allele1="X"
if(fil_format(1:8)=="ldmipgen")then
   read(4,*) !skip first line
   if(ifilnam2==0)then
     ioption=index(file_name,'.gen')
     file_name(ioption:ioption+1)='.r'  !.ren file
     open(5,file=file_name,status='unknown')
     print *,' renumbered genotypes: ',File_name
   else
     line2="head -1 "; line2(9:)=trim(genotypefile); line2(len_trim(line2)+1:)=" >w.w"
     call system(trim(line2),istatus)
     if(istatus/=0)print *,' ERR: ',line2
     line2="gawk 'BEGIN{getline <#w.w#;for(i=2;i<=NF;i++)lst[$i]=(i-1)};NR>1{print $1,lst[$2],$2}' "  !filnam2 > gghatvr3.tmp"
     line2(22:26)='"w.w"';  line2(len_trim(line2)+1:)=trim(mapfile); line2(len_trim(line2)+1:)=" > gghatvr3.tmp"
     call system(trim(line2),istatus)
     if(istatus/=0)print *,' ERR: ',line2
     open(9,file="gghatvr3.tmp",status='old')
     nchr=0; nloc_chr=0
     do i=1,10000000
       read(9,*,end=8)ichr,iloc
       nloc_chr(ichr)=nloc_chr(ichr)+1
       iorder(i)=iloc
       nchr=max(ichr,nchr)
     enddo
8    close(9)
    do i=1,nchr
      line2=trim(genotypefile)
      if(i<10)then
        write(line2(len_trim(line2)+1:len_trim(line2)+1),'(i1)')i  
      else
        write(line2(len_trim(line2)+1:len_trim(line2)+2),'(i2)')i  
      endif
      open(10+i,file=trim(line2),status='unknown')
     enddo
   endif
endif
!  read(4,*,end=9)ids(ncount),(iifix(1,j,ncount),j=1,nloc) !imean is actually ID number  

nhap=0; ifix=0; wt=1.; imean=1; no_means=1; frq=0.; nfrq=0;ncount=0; 
do 
 ncount=ncount+1
 if(fil_format(1:8)=="ldmipgen" .or. fil_format(1:8)=="ldmipmrk" .or. fil_format(1:3)=="012" .or. fil_format(1:3)=="bgl" .or. fil_format(1:5)=="plink")then !OPTION0
  if(fil_format(1:8)=="ldmipgen")then   !OPTION1
   read(4,'(a)',end=9)line
   call readline(trim(line),nloc,ids(ncount),iifix(1:2,1:nloc,ncount),allele1)
   if(ifilnam2==0)then
     write(5,'(i10,1000000i2)')ids(ncount),(iifix(1,j,ncount),j=1,nloc) !imean is actually ID number
     write(5,'(i10,1000000i2)')ids(ncount),(iifix(2,j,ncount),j=1,nloc)   
   else
     ll=0
     do ichr=1,nchr
       write(10+ichr,'(i10,1000000i2)')ids(ncount),(iifix(1,iorder(j),ncount),j=ll+1,ll+nloc_chr(ichr))
       write(10+ichr,'(i10,1000000i2)')ids(ncount),(iifix(2,iorder(j),ncount),j=ll+1,ll+nloc_chr(ichr))
       ll=ll+nloc_chr(ichr)
     enddo
   endif
  elseif(fil_format(1:3)=="bgl")then  !OPTION1
   if(ilocus==0)then
      file_name=trim(file_name)//".ren"
      close(55,iostat=ios) !close(55) if open       
      open(55,file=file_name,status='unknown')
      print *,' renumbered genotypes: ',File_name
   endif
   read(4,'(a)',end=9)line
   call readline2(trim(line),nanim,ilocus,iid,iifix(1:2,1:nloc,1:pndat),allele1)
   write(55,'(i10,1000000i2)')iid,((iifix(i,ilocus,j),i=1,2),j=1,pndat) !imean is actually ID number
  elseif(fil_format(1:3)=="012")then !OPTION1
     read(4,*,end=9)ids(ncount),(iifix(1,j,ncount),j=1,nloc)
  elseif(fil_format(1:5)=="plink")then
     read(4,end=9)ifour(1:nbytes)
     ianim=0; sumd=0.0
     do j=1,nbytes
        jval=ibits(ifour(j),jpos,2)
        do k=1,4
           IZt=0; ianim=ianim+1
           if(ianim>pndat)exit
           afix(ncount,ianim)=2.  !2 is later substracted
           if(jval(k)==1)then      !note homozygot has automatically value '0'
              print *,' ERROR: no missing genotypes allowed; missing genotypes are aribitrarily set to heterozygot ' ; IZt=1
           elseif(jval(k)==2)then  !change heterozygot to '1'
              afix(ncount,ianim)=1+2.
           elseif(jval(k)==3)then  !change homozygot to '2'
              afix(ncount,ianim)=2+2.
           endif
           frq(ncount)=frq(ncount)+afix(ncount,ianim)-2.
           nfrq(ncount)=nfrq(ncount)+2
        enddo
     enddo
  else                           !OPTION1
   read(4,*,end=9)ids(ncount),(iifix(1,j,ncount),j=1,nloc) !imean is actually ID number
   read(4,*)ids(ncount),(iifix(2,j,ncount),j=1,nloc)
  endif                          !OPTION1_END
  icum=0
  do j=1,nloc
   if(fil_format(1:3)=="012")then  !0/1/2 genotype coding
    if(iifix(1,j,ncount)>=0 .and. iifix(1,j,ncount)<=2)then
      afix(j,ncount)=iifix(1,j,ncount)+2.   !2 is later substracted
      frq(j)=frq(j)+afix(j,ncount)-2.
      nfrq(j)=nfrq(j)+2
     else
       afix(j,ncount)=0
     endif
    endif
  enddo

 elseif(fil_format(1:8)=="ldmipout")then   !OPTION0
  read(4,*,end=9)ianim,iloc,gen1,gen2,gen3,gen4 
     if(ianim>pndat)cycle
     if(iloc>nloc)cycle
    ifix(1,iloc,ianim)=gen3+gen4
    ifix(2,iloc,ianim)=gen2+gen4
    frq(iloc)=frq(iloc)+gen2+gen3+2.*gen4
    nfrq(iloc)=nfrq(iloc)+2
 endif                                !OPTION0_END
end do                                !END_LOOP_READ(4)
9 continue
close(4); 
if(fil_format(1:3)/="bgl")close(5,iostat=ios) 
nhap=2 !always
if(fil_format(1:8)=="ldmipgen" .or. fil_format(1:8)=="ldmipmrk" .or. fil_format(1:3)=="bgl")then !1/2 allele coding
  do j=1,nloc
    do i=1,pndat
     if(iifix(1,j,i)>0 .and. iifix(1,j,i)<=2)then
      afix(j,i)=iifix(1,j,i)+iifix(2,j,i)
      frq(j)=frq(j)+afix(j,i)-2.
      nfrq(j)=nfrq(j)+2
     else
       afix(j,ncount)=0
     endif
    enddo
 enddo
endif



l=nloc
if(calc_freqs)then
  where(nfrq(1:nloc)>0)frq(1:nloc)=frq(1:nloc)/nfrq(1:nloc)
else
  do i=1,nloc
    read(80,*)frq(i)
  enddo
endif
where(frq==.5)frq=.499   !because afix=0.0 is used to denote missing records and frq=.5 can result afix=0.0
print *,'frq(1:5)',frq(1:5)

if(fil_format(1:8)=="ldmipout")then
   denom=sum(frq(1:l)*(1.-frq(1:l)))
   if(abs(ivanraden)==1)then
     forall(ii=1:l,jj=1:pndat)ifix(:,ii,jj)=(ifix(:,ii,jj)-frq(ii))  !/sqrt(frq(ii)*(1.-frq(ii)))
   else
     forall(ii=1:l,jj=1:pndat)ifix(:,ii,jj)=(ifix(:,ii,jj)-frq(ii))/sqrt(frq(ii)*(1.-frq(ii)))
   endif
else !ldmipmrk
   denom=2.*sum(frq(1:l)*(1.-frq(1:l)))
   print *,'afix',afix(1:5,1)
   print *,'afix',afix(1,1:5)
   do i=1,pndat
    if(abs(ivanraden)==1)then
      where(afix(1:l,i)>0 .and. frq(1:l)>0 .and. frq(1:l)<1.)
         afix(1:l,i)=(afix(1:l,i)-2.-2.*frq(1:l))  !/sqrt(2.*frq(1:l)*(1.-frq(1:l)))
      elsewhere
         afix(1:l,i)=0.
      endwhere
    else
      where(afix(1:l,i)>0 .and. frq(1:l)>0 .and. frq(1:l)<1.)
         afix(1:l,i)=(afix(1:l,i)-2.-2.*frq(1:l))/sqrt(2.*frq(1:l)*(1.-frq(1:l)))
      elsewhere
         afix(1:l,i)=0.
      endwhere
    endif
   enddo
   print *,'afix2',afix(1:5,pndat)
   print *,'afix2',afix(1,1:5)
endif
if(abs(ivanraden)==1)then
  denomsum=denomsum+denom
else
  denomsum=denomsum+nloc
endif
print *,' No of loci,denom     =',nloc,denom

if(igiv>=0)then
l=nloc
!!$OMP PARALLEL DO PRIVATE(i,j,ii,jj) SHARED(pndat,l,gmat,ifix,afix,option,ivanraden,frq)
do i=1,pndat  !,1,-1 !pndat  !reverse order so that last row takes little time
 ii=pndat+i
 do j=1,i-1
 if(fil_format(1:8)=="ldmipout")then
   jj=pndat+j
   Gmat(i,j)=Gmat(i,j)+dot(ifix(1,1:l,i),ifix(1,1:l,j))
   gmat(i,jj)=Gmat(i,jj)+dot(ifix(1,1:l,i),ifix(2,1:l,j))
   gmat(ii,j)=Gmat(ii,j)+dot(ifix(2,1:l,i),ifix(1,1:l,j))
   gmat(ii,jj)=Gmat(ii,jj)+dot(ifix(2,1:l,i),ifix(2,1:l,j))
 else
   Gmat(j,i)=Gmat(j,i)+dot(afix(1:l,i),afix(1:l,j))
   if(abs(ivanraden)==1)then
     Gmat(pndat+j,pndat+i)=Gmat(pndat+j,pndat+i)+2.*sum(frq(1:l)*(1.-frq(1:l)),mask=(afix(1:l,i)/=0.0 .and. afix(1:l,j)/=0.0))
   else
     Gmat(pndat+j,pndat+i)=Gmat(pndat+j,pndat+i)+count(afix(1:l,i)/=0.0 .and. afix(1:l,j)/=0.0)
   endif
 endif !option
 enddo
 if(fil_format(1:8)=="ldmipout")then
   Gmat(i,i)=Gmat(i,i)+dot(ifix(1,1:l,i),ifix(1,1:l,i))
   Gmat(ii,ii)=Gmat(ii,ii)+dot(ifix(2,1:l,i),ifix(2,1:l,i))
   Gmat(i,ii)=Gmat(i,ii)+dot(ifix(1,1:l,i),ifix(2,1:l,i))
   Gmat(ii,i)=Gmat(i,ii)
 else
   Gmat(i,i)=Gmat(i,i)+dot(afix(1:l,i),afix(1:l,i))
   if(abs(ivanraden)==1)then
     Gmat(pndat+i,pndat+i)=Gmat(pndat+i,pndat+i)+2.*sum(frq(1:l)*(1.-frq(1:l)),mask=(afix(1:l,i)/=0.0))
   else
     Gmat(pndat+i,pndat+i)=Gmat(pndat+i,pndat+i)+count(afix(1:l,i)/=0.0)
   endif
   if(i==pndat)print *,'Gmat(n,n)',i,Gmat(i,i),Gmat(pndat+i,pndat+i)
 endif
enddo
!!$OMP END PARALLEL DO
endif
!Gmat(1:pndat,1:pndat)=Gmat(1:pndat,1:pndat)+matmul(transpose(afix(1:l,1:pndat)),afix(1:l,1:pndat))
nnloc=nnloc+nloc

!if(igiv==-2)then  !MOVE TO WITHIN CHROM AREA
!!  print *,' start cluster analysis'
!!  if(option(1:8)=="ldmipout")afix(1:L,1:pndat)=ifix(1,1:L,1:pndat)+ifix(2,1:L,1:pndat)
!!  call clust(afix,r2crit,ichrom)
!endif
ENDDO !ICHROM
if(igiv<0)stop' Finished (without setting up G)'

if(fil_format(1:8)=="ldmipout")then
  Gmat=Gmat/denomsum  !divide by sum of heterozygosities
else
  where(Gmat(pndat+1:2*pndat,pndat+1:2*pndat)>0.0)Gmat(1:pndat,1:pndat)=Gmat(1:pndat,1:pndat)/Gmat(pndat+1:2*pndat,pndat+1:2*pndat)
  print *,"Gmat",Gmat(1:5,1)
endif

IF(fil_format(1:8)=="ldmipout")then  !account for uncertainties in genotype probs
wt=1.               !scale back the too high estimates of the diagonal
forall(ii=1:pndat*2)wt(ii)=Gmat(ii,ii)
wt2=1.-wt
wt=sqrt(wt)
where(wt<1.)wt=1.
where(wt2<=0.)wt2=0.0
wt2=sqrt(wt2)
forall(ii=1:pndat*2,jj=1:pndat*2)
    Gmat(ii,jj)=Gmat(ii,jj)/(wt(ii)*wt(jj))
    Amat(ii,jj)=Amat(ii,jj)*wt2(ii)*wt2(jj)
endforall
!forall(ii=1:pndat*2,wt(ii),ii)<1.)Gmat(ii,ii)=1.   !set the diagonals that are too small to 1
Gmat=Gmat+Amat
ENDIF

open(3,file='gghatvr3.g',status='unknown')
do i=1,pndat
ii=i+pndat
do j=1,i
  jj=j+pndat
 if(fil_format(1:8)=="ldmipout")then
    gmat(i,j)=(Gmat(i,j)+Gmat(ii,j)+Gmat(i,jj)+Gmat(ii,jj))/2. 
    gmat(j,i)=gmat(i,j)
    write(3,*)i,j,gmat(i,j)
 else
    gmat(i,j)=gmat(j,i)
    write(3,*)i,j,gmat(i,j),ids(i),ids(j)   
 endif
enddo
enddo
close(3)
print *,' total number of loci',nnloc
print *," Gmatrix written to 'gghatvr3.g'"



!inversion
23 continue
IF(IGIV>0)THEN
  if(allocated(gmat))deallocate(gmat)                     !reduce memory usage
  if(allocated(amat))deallocate(amat)  
  allocate(gmat(pndat,pndat))
  open(3,file='gghatvr3.g',status='old')
  do i=1,pndat
     do j=1,i
        read(3,*)ii,jj,gmat(i,j)
        gmat(j,i)=gmat(i,j)
     enddo
     gmat(i,i)=gmat(i,i)+add_diag(1)   !stabilise matrix
  enddo
  close(3)
  !   forall(ii=1:pndat)Gmat(ii,ii)=Gmat(ii,ii)+add_diag(1)  !stabilise matrix
  
CALL POTRF(gmat(1:pndat,1:pndat),info=inf) !"L", n, a, ia, ja, desc_a, info)
if(inf/=0)then
  print *,' ERROR: POTRFinfo=',inf
  stop
endif
det_log=0.0
do i=1,pndat
  det_log=det_log+2.*log(gmat(i,i))
enddo
print *,' log(det(G)) =  ',det_log
CALL POTRI(gmat(1:pndat,1:pndat),info=inf) 
if(inf/=0)print *,' POTRIinfo=',inf

      open(3,file='gghatvr3.giv',status='unknown')
      DO I=1,PNDAT;
      do j=1,i
         aa=gmat(j,i)
         if(fil_format(1:8)=="ldmipout")then
           WRITE(3,*)i,j,aa
         else
           WRITE(3,*)i,j,aa,ids(i),ids(j)
         endif
         gmat(i,j)=aa
      ENDDO
      enddo
print *," inverse written to 'gghatvr3.giv'"
ENDIF !IGIV



!Gmat(1:pndat,1:pndat)=inv(Gmat(1:pndat,1:pndat))                    
!
!open(3,file='ghat.giv',status='unknown')
!do i=1,pndat
!do j=1,i
!  write(3,*)i,j,Gmat(i,j)
!enddo
!enddo
!close(3)
!print *," Ginvers written to 'ghat.giv'"

contains
subroutine  readline(line,nloc,iid,iifix,allele1)
character(len=*) :: line
character(len=1) :: allele1(:)
integer :: iloc,lt,nloc,iid,iifix(:,:),i,j,iread,imat

   iloc=0; LT=len_trim(line)
   i=1
   DO IREAD=0,NLOC*2
     j=verify(line(i:)," ")
    if(j>0)i=i+j-1
    if(iread==0)then
      read(line(i:),*)iid;
    else
     imat=2
     if(mod(iread,2)==1)then
       iloc=iloc+1; imat=1
     endif
     if(line(i:i)/="0")then  !nonmissing
       if(allele1(iloc)=="X")allele1(iloc)=line(i:i)  !decide allele1 is coded as 1; any other allele will be coded 2
       if(line(i:i)==allele1(iloc))then
         iifix(imat,iloc)=1
       else
         iifix(imat,iloc)=2
       endif
     else
         iifix(imat,iloc)=0
     endif
   endif
   j=max(scan(line(i:),' '),2) !move at least 1
   I=i+j-1
   if(i>LT)exit
  ENDDO
if(iloc/=nloc)print *,' ERR read ',iloc,' loci' 
end subroutine

subroutine  readline2(line,nanim,ilocus,iid,iifix,allele1)
character(len=*) :: line
character(len=1) :: allele1(:)
integer :: ilocus,iloc,lt,nloc,iid,iifix(:,:,:),i,j,iread,imat,nanim
ilocus=ilocus+1
nanim=size(iifix,3)

   iloc=0; LT=len_trim(line); !iloc is actually animal
   i=2  !skip 'M'
   DO IREAD=0,NANIM*2
     j=verify(line(i:)," ")
    if(j>0)i=i+j-1
    if(iread==0)then
      read(line(i:),*)iid;
    else
     imat=2
     if(mod(iread,2)==1)then
       iloc=iloc+1; imat=1
     endif
     if(line(i:i)/="0")then  !nonmissing
!       if(iread<20)print '(2a)','allele1 ',allele1(ilocus)
       if(allele1(ilocus)=="X")allele1(ilocus)=line(i:i)  !decide allele1 is coded as 1; any other allele will be coded 2
       if(line(i:i)==allele1(ilocus))then
         iifix(imat,ilocus,iloc)=1
       else
         iifix(imat,ilocus,iloc)=2
       endif
     else
         iifix(imat,ilocus,iloc)=0
     endif
   endif
   j=max(scan(line(i:),' '),2) !move at least 1
   I=i+j-1
   if(i>LT)exit
   ENDDO
if(iloc/=nanim)print *,' ERR readline2: ',iloc,' animals read' 
!print *,line(1:80)
!print *,'line2,',iifix(1,ilocus,1:10)
!print *,'line2,',iifix(2,ilocus,1:10)
end subroutine


end program
