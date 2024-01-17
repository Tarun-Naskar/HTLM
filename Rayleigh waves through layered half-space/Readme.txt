# HTLM_R

Manuscript title –Efficient Modeling of Surface Wave Dispersion with Higher-Order Thin Layer Method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Program:: HTLM_R.m

 Coded by: Mrinal Bhaumik and Tarun naskar
 Indian Institute of Technology Madras, India

 Last revision date:
 16/01/2024

Inform issues to:  wavepropagationlab@gmail.com


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SCRIPTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

"HTLM_R.m"


Description: 

"HTLM_R.m" is a Matlab script to generate dispersion curves of Rayleigh wave using higher order thin layer method.
Users are encouraged  to refer to the corresponding research paper for a detailed explanation of the underlying concepts and methods.

The code HTLM implements higher order elements, which makes the traditional TLM faster and more accurate.

Below are descriptions for the scripts, the functions, and the example.
__________________________________________________________________________
__________________________________________________________________________




To run the code:

Step 1 - Prepare an Excel sheet containing plate profile information as follows:

       Data Format
       h | vs | vp/mu | rho

       where h - layer thickness in m, vs - shear wave velocity in m/s, vp/mu - compression wave velocity in m/s or Poisson’s ratio, 
       and rho – density in kg/m^3.

       Four examples of Excel sheet is provided with the code.


Step 2 – Assign the parameters. 

 	fmax     – maximum frequency (Hz) up to which the dispersion spectra are required.
 	df       – the frequency resolution (Hz).
 	d        – is the order of HTLM elements. The method can take from linear (d=1) to 15th order (d = 15).
 	dh       – is the thickness (m) of the thin discretized element
 	PML      – Number of PMDL layer require to model the half-space. For soil model (m scale)
                   10 PMDL layer is sufficient. For crustal model (km scale) 15 is recommended.
 	fm       – is the frequency (Hz) where mode shape is required. Other wise keep it empty [].
 	mode_n   - mode index.

__________________________________________________________________________

*** (d = 1, corresponding to existing linear element TLM. dh should be at 
least L_min/8 to L_min / 10. For higher order elements dh can be selected 
between L_min/3-L_min/4).


All the necessary calculations are performed in the “Main_code” sub-file.

__________________________________________________________________________

Output:
	v        – is the matrix containing phase velocities of all the possible modes
 	ev_mat   – is the matrix containing mode shapes at frequency “fm”
 	dof      – provide the total degree of freedom (except half-space), involves in the calculation.






%%%%%%%%%%%%%%%%%%%%%%%%            EXAMPLES           %%%%%%%%%%%%%%%%%%%%%%%%%%

Example: 1


read'Profile_I.xlsx'
 
fmin    = 1;                    % Minimum required frequency
fmax    = 50;                   % Maximum required frequency
df      = 0.5;                  % Frequency resolution
d       = 6;                    % Order of HTLM (1-15),(1-Linear, 2-Quadratic, 3-Cubic, 4-Quartic)
dh      = 10;                   % Thickness of the thin layers.

PML     = [10] ;                % No of PMDL layer (default value 10)

fm      = [30];                 % Plot mode shape at frequency "fm"; else, fm = [], fm <= fmax
mode_n  = [3];                  % mode index. plots nth mode at frequency "fm"

w       = fmin : df : fmax;     % Frequency scale


Output:

 Figure 1 (Dispersion curves)
 Figure 2 (3rd mode shape of 30 Hz)

