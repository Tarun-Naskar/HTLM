# HTLM_L

Manuscript title –Higher-order thin layer method as an efficient forward model for calculating 
dispersion curves of surface and Lamb waves in layered media

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Program:: HTLM_L.m

 Coded by: Mrinal Bhaumik and Tarun naskar
 Indian Institute of Technology Madras, India

 Last revision date:
 16/01/2024

Inform issues to:  wavepropagationlab@gmail.com


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SCRIPTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

"HTLM_L.m"


Description: 

"HTLM_L.m" is a Matlab script to generate dispersion spectra of Love wave (SH-wave or antiplane shear) using higher order thin layer method.
Users are encouraged  to refer to the corresponding research paper for a detailed explanation of the underlying concepts and methods.

The code HTLM implements higher order elements, which makes the traditional TLM faster and more accurate.

Below are descriptions for the scripts, the functions, and the examples.
__________________________________________________________________________
__________________________________________________________________________




To run the code:

Step 1 - Prepare an Excel sheet containing plate profile information as follows:

       Data Format
       h | vs | rho

       where h - layer thickness in m, vs - shear wave velocity in m/s, 
       and rho – density in kg/m^3.

       Four examples of Excel sheet is provided with the code.


Step 2 – Assign the parameters. 

    fmin     – minimum frequency (Hz).
    fmax     – maximum frequency (Hz) up to which the dispersion spectra are required.
    df       – the frequency resolution (Hz).
    d        – is the order of HTLM elements. The method can take from linear (d=1) to 15th order (d = 15).
    dh       – is the thickness (m) of the thin discretized element
    PML      – Number of PMDL layer require to model the half-space. 
               For soil model (m scale) 10 PMDL layer is sufficient. For crustal model (km scale) 15 is recommended.

__________________________________________________________________________

*** (d = 1, corresponding to existing linear element TLM. dh should be at 
least L_min/8 to L_min / 10. For higher order elements dh can be selected 
between L_min/3-L_min/4).


All the necessary calculations are performed in the “Main_code” sub-file.

__________________________________________________________________________

Output:
   v        – is the matrix containing phase velocities of all the possible modes.
   EV       - eigenvectors or mode shape (3D matrix : depth * mode no * frequency ) 
   dh_vec   – depth of the nodes.






%%%%%%%%%%%%%%%%%%%%%%%%            EXAMPLES           %%%%%%%%%%%%%%%%%%%%%%%%%%

Example: 1


read'Profile_II.xlsx'
 
fmin    = 0;                    % minimum frequency
fmax    = 50;                   % Maximum required frequency
df      = 0.5;                  % Frequency resolution

d       = 6;                    % Order of HTLM (1-15),(1-Linear, 2-Quadratic, 3-Cubic, 4-Quartic)
dh      = 10;                   % Thickness of the thin layers.

PML     = 10 ;                  % No of PMDL layer (default value 10)
w       = fmin : df : fmax;     % Frequency scale


Output:

 Figure 1 (Dispersion spectra)

