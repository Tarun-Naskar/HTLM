# HTLM_Lamb

Manuscript title – Higher-order thin layer method as an efficient forward model for calculating 
dispersion curves of surface and Lamb waves in layered media


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Program:: HTLM_plate.m

 Coded by: Mrinal Bhaumik and Tarun naskar
 Indian Institute of Technology Madras, India

 Last revision date:
 16/01/2024

Inform issues to:  wavepropagationlab@gmail.com


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SCRIPTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

"HTLM_Lamb.m"


Description: 

"HTLM_Lamb.m" is a Matlab script to generate dispersion curves of Lamb wave using higher order thin layer method. Users are encouraged 
to refer to the corresponding research paper for a detailed explanation of the underlying concepts and methods.

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

       An example of Excel sheet is provided with the code.


Step 2 – Assign the parameters. 

    fmin     – Starting frequency (Hz).
    fmax     – maximum frequency (Hz) up to which the dispersion spectra are required.
    df       – the frequency resolution (Hz).
    d        – is the order of HTLM elements. The method can take from linear (d=1) to 15th order (d = 15).
    dh       – is the thickness (m) of the thin discretized element
    fm       – is the frequency (Hz) where mode shape is required. Otherwise keep it empty [].
    mode_n   - mode index.

__________________________________________________________________________

*** (d = 1, corresponding to existing linear element TLM. dh should be at 
least L_min/8 to L_min / 10. For higher order elements dh can be selected 
between L_min/3-L_min/4).


All the necessary calculations are performed in the “Main_code” sub-file.

__________________________________________________________________________

Output:
     V      – is the matrix containing phase velocities of all the possible modes
     ev_mat – is the matrix containing mode shapes at frequency “fm”
     dof    – provide the total degree of freedom (except half-space), involved in the calculation.






%%%%%%%%%%%%%%%%%%%%%%%%            EXAMPLES           %%%%%%%%%%%%%%%%%%%%%%%%%%

Example: 1


read'Profile_Plate.xlsx'
 
fmin    = 0;                    % Starting frequency
fmax    = 30000;                % Maximum required frequency
df      = 100;                  % Frequency resolution
d       = 6;                    % Order of HTLM (1-15),(1-Linear, 2-Quadratic, 3-Cubic, 4-Quartic)
dh      = 0.05;                 % Thickness of the thin layers. 
w       = fmin : df : fmax;     % Frequency scale
fm      = [20000];              % Plot mode shape at frequency "fm"; else, fm = [], fm <= fmax
mode_n  = [2];                  % mode index. plots nth mode at frequency "fm"


Output:

 Figure 1 (Dispersion spectra)
 Figure 2 (2nd mode shape of 20000 Hz)
