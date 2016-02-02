C-Stores the values of the two types of homogeneous
C-solutions required to construct the global solution
C-within a given subdomain. First index represents
C-type of subdomain. There are as many types of subdomains
C-as there are different values of subdomain heights.
      common /homogsol/slhzetauvl(nsdtype,nxhp,ny,nzloc),
     >slhzetauvr(nsdtype,nxhp,ny,nzloc),
     >slhwl(nsdtype,nxhp,ny,nzloc),
     >slhwr(nsdtype,nxhp,ny,nzloc),
     >slhtempl(nsdtype,nxhp,ny,nzloc),
     >slhtempr(nsdtype,nxhp,ny,nzloc)
