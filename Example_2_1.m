%% Exercise 2.1
clear all 

%% Parameters
m=1;                                % mass
k=10000;                            % stiffness
wn=sqrt(k/m);                       % undamped natural frequency
z=0.1;                              % damping ratio (set to 0.001, 0.01, 0.1, 0.999)
c=2*z*wn*m;                         % damping coefficient
wd=sqrt(1-z^2)*wn;                  % damped natural frequency

%% Time Vector
dt=0.001;                           % time resolution
T=100;                              % duration of time signal
t=0:dt:T;                           % time vector

%% Impulse response
h=1/(m*wd)*exp(-z*wn*t).*sin(wd*t); % displacement IRF
 
%% Plot results
plot(t,h,'linewidth',2,'Color',[1 1 1]*0.2);    
axis([0,0.8,-0.01,0.01]);
grid;axis square
set(gca,'fontsize',16)
xlabel('time (s)');
ylabel('displacement IRF (m/Ns)');

%% function

/POST1
mode_data=
matid = 
allsel
ESEL,S,ENAME,,214
*get, num_els, ELEM, 0, COUNT
*get, el_min, ELEM, 0, NUM, MIN
*get,nb_modes,active,0,set, nset

nbmat=85

nb_modes = nb_modes/4

*DIM,mode_data,ARRAY,nb_modes,19  
*dim,matid,string,80,nbmat  

*do,id,1,nbmat
matid(1,id) = strcat(strcat('mat',chrval(id)), ' ')
*enddo

*stat,matid

allsel

*get, mass_x, ELEM,,MTOT,X 
*get, mass_y, ELEM,,MTOT,Y 
*get, mass_z, ELEM,,MTOT,Z 
*stat, mass_x
*stat, mass_y
*stat, mass_z

*DO,i,1,nb_modes
	mode_data(i,1) = i
    set,1,2*i-1,,1
	*get,freq,mode,2*i-1,dfrq
	mode_data(i,2) = freq
	el=el_min
	allsel
	PRENERGY, SENE,
	*get, tenem, PRENERGY,, TOTE, 1
	ESEL,S,ENAME,,214
	ETABLE,SENE_BEARING,SENE
	*GET,se_brg1,ELEM,el,ETAB,SENE_BEARING
	mode_data(i,3) = se_brg1/tenem*100
	el = elnext(el)
	*GET,se_brg2,ELEM,el,ETAB,SENE_BEARING
	mode_data(i,4) = se_brg2/tenem*100
	allsel

	ESEL,S,MAT,,mat1
	*do,j,2,nbmat
	    ESEL,A,MAT,,%matid(1,j)%
	*enddo
	etab, sene_rotor, SENE
	sabs,1
	ssum
	*GET,sesum_rotor,SSUM,,ITEM,sene_rotor
	*stat,sesum_rotor
	mode_data(i,5) = sesum_rotor/tenem*100
	allsel
	
	ESEL,S,MAT,,mat1
	*do,j,2,18
	    ESEL,A,MAT,,%matid(1,j)%
	*enddo
	ESEL,A,MAT,,mat83
	etab, sene_front_shaft, SENE
	sabs,1
	ssum
	*GET,sesum_front_shaft,SSUM,,ITEM,sene_front_shaft
	*stat,sesum_front_shaft
	mode_data(i,6) = sesum_front_shaft/tenem*100
	
	
	allsel	
	ESEL,S,MAT,,mat84
	etab, sene_hpc1, SENE
	sabs,1
	ssum
	*GET,sesum_hpc1,SSUM,,ITEM,sene_hpc1
	*stat,sesum_hpc1
	mode_data(i,7) = sesum_hpc1/tenem*100
	
	
	allsel	
	ESEL,S,MAT,,mat85
	etab, sene_hpc2, SENE
	sabs,1
	ssum
	*GET,sesum_hpc2,SSUM,,ITEM,sene_hpc2
	*stat,sesum_hpc2
	mode_data(i,8) = sesum_hpc2/tenem*100
	
	allsel
	ESEL,S,MAT,,mat19
	*do,j,20,25
	    ESEL,A,MAT,,%matid(1,j)%
	*enddo
	etab, sene_hpc3_4, SENE
	sabs,1
	ssum
	*GET,sesum_hpc3_4,SSUM,,ITEM,sene_hpc3_4
	*stat,sesum_hpc3_4
	mode_data(i,9) = sesum_hpc3_4/tenem*100	
	
	allsel	
	ESEL,S,MAT,,mat26
	ESEL,A,MAT,,mat27
	etab, sene_hpc5, SENE
	sabs,1
	ssum
	*GET,sesum_hpc5,SSUM,,ITEM,sene_hpc5
	*stat,sesum_hpc5
	mode_data(i,10) = sesum_hpc5/tenem*100

	allsel
	ESEL,S,MAT,,mat28
	*do,j,29,35
	    ESEL,A,MAT,,%matid(1,j)%
	*enddo
	etab, sene_impeller, SENE
	sabs,1
	ssum
	*GET,sesum_impeller,SSUM,,ITEM,sene_impeller
	*stat,sesum_impeller
	mode_data(i,11) = sesum_impeller/tenem*100	

	allsel	
	ESEL,S,MAT,,mat36
	ESEL,A,MAT,,mat37
	etab, sene_preloadshaft, SENE
	sabs,1
	ssum
	*GET,sesum_preloadshaft,SSUM,,ITEM,sene_preloadshaft
	*stat,sesum_preloadshaft
	mode_data(i,12) = sesum_preloadshaft/tenem*100		

	allsel
	ESEL,S,MAT,,mat38
	*do,j,39,54
	    ESEL,A,MAT,,%matid(1,j)%
	*enddo
	etab, sene_midstubshaft, SENE
	sabs,1
	ssum
	*GET,sesum_midstubshaft,SSUM,,ITEM,sene_midstubshaft
	*stat,sesum_midstubshaft
	mode_data(i,13) = sesum_midstubshaft/tenem*100		
	
	allsel
	ESEL,S,MAT,,mat55
	*do,j,56,68
	    ESEL,A,MAT,,%matid(1,j)%
	*enddo
	etab, sene_hpt1, SENE
	sabs,1
	ssum
	*GET,sesum_hpt1,SSUM,,ITEM,sene_hpt1
	*stat,sesum_hpt1
	mode_data(i,14) = sesum_hpt1/tenem*100	

	allsel
	ESEL,S,MAT,,mat69
	*do,j,70,73
	    ESEL,A,MAT,,%matid(1,j)%
	*enddo
	etab, sene_hpt2, SENE
	sabs,1
	ssum
	*GET,sesum_hpt2,SSUM,,ITEM,sene_hpt2
	*stat,sesum_hpt2
	mode_data(i,15) = sesum_hpt2/tenem*100	
	
	
	allsel
	ESEL,S,MAT,,mat74
	*do,j,75,80
	    ESEL,A,MAT,,%matid(1,j)%
	*enddo
	etab, sene_rearnut, SENE
	sabs,1
	ssum
	*GET,sesum_rearnut,SSUM,,ITEM,sene_rearnut
	*stat,sesum_rearnut
	mode_data(i,16) = sesum_rearnut/tenem*100	

	allsel	
	ESEL,S,MAT,,mat81
	ESEL,A,MAT,,mat82
	etab, sene_tieshaft, SENE
	sabs,1
	ssum
	*GET,sesum_tieshaft,SSUM,,ITEM,sene_tieshaft
	*stat,sesum_tieshaft
	mode_data(i,17) = sesum_tieshaft/tenem*100		
	
	allsel
	
	*GET,part_fact_z,MODE,i,PFAC,,DIREC,z
    *stat, part_fact_z
    effm_tz = part_fact_z*part_fact_z
    *stat, effm_tz
    mode_data(i,18)= effm_tz/mass_z*100
	
	allsel
	
	*GET,part_fact_y,MODE,i,PFAC,,DIREC,y
    *stat, part_fact_y
    effm_ty = part_fact_y*part_fact_y
    *stat, effm_ty
    mode_data(i,19)= effm_ty/mass_y*100
	
	*stat,i
	
*ENDDO
*STAT,mode_data

*mwrite,mode_data,modes,csv,,jik
(20(f15.6,','))


*CFOPEN,modes_detailed,txt
*VWRITE,'Mode Num','Freq(Hz)','Bea1 SE(%)','Bea2 SE(%)','Rotor SE(%)','Front Shaft SE(%)','HPC1 SE(%)','HPC2 SE(%)','HPC3-4 SE(%)','HPC5 SE(%)','Impeller SE(%)','PreloadS SE(%)','MidstubS SE(%)','HPT1 SE(%)','HPT2 SE(%)','Rearnut SE(%)','Tieshaft SE(%)','Eff Mass Z (%)','Eff Mass Y (%)'
%16C%16C%16C%16C%16C%16C%16C%16C%16C%16C%16C%16C%16C%16C%16C%16C%16C%16C%16C

*VWRITE,mode_data(1,1),mode_data(1,2),mode_data(1,3),mode_data(1,4),mode_data(1,5),mode_data(1,6),mode_data(1,7),mode_data(1,8),mode_data(1,9),mode_data(1,10),mode_data(1,11),mode_data(1,12),mode_data(1,13),mode_data(1,14),mode_data(1,15),mode_data(1,16),mode_data(1,17),mode_data(1,18),mode_data(1,19)                       
%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G%16.8G
*CFCLO
