% This is a script that uses Lead-DBS's "ea_map_coords" to transform estimated microelectrode
% locations from native space (anat_t1.nii) to MNI space (glanat.nii).

%% For this script to work, the following steps have to precede:
% 1) pre-process MRI data in Lead-DBS
% 2) run "x2_trajcalc.R" to estimated microelectrode locations in native space

%% loop through exploratory and permanenet electrode coordinates
%% and transform each coordinate to MNI space
dirs = { '_nogithub/coords/coords_expl.csv' '_nogithub/coords/coords_check.csv' };

% the loop itself
for x = 1:2
    nat2mni( dirs{x}, dirs{x} )
end


%% prepare a function that:
% 1) takes-in a pre-processed coordinates file (input_file),
% 2) computes raw2mni transformation of input coordinates,
% 3) saves the result as output_file (might be the same as input_file)
function d = nat2mni( input_file, output_file )

    % read the table with estimated coordinates in native space
    d = readtable( input_file );
    nrow = height(d); % extract the number of rows
    d( :, {'mni_x','mni_y','mni_z'} ) = array2table( zeros(nrow, 3) ); % prepare new columns
    
    % run a loop through all rows of d
    for i = 1:nrow
    
        % set-up the directory to patient's folder
        dir = strcat( '_nogithub/mri/lead/', convertCharsToStrings( d.id(i) ) );
        
        % prepare the native coordinates in voxels
        c = table2array( d( i, {'nat_x','nat_y','nat_z'} ) ); % the coordinates
        n = char( strcat( dir, '/anat_t1.nii' ) ); % native space
    
        % skip missing contacts or non-ANTs transformations (for now)
        if sum( isnan(c) ) > 0 || d.norm_algo(i) == "FSL_FNIRT"
            d( i, {'mni_x','mni_y','mni_z'} ) = array2table( [ NaN NaN NaN ] );
        else
            c_vox = ea_mm2vox( c , n ); % mapping from mm to voxels in native space
            t = char( strcat( dir, '/y_ea_inv_normparams.nii' ) ); % normalisation inverse transform
            d( i, {'mni_x','mni_y','mni_z'} ) = array2table( ea_map_coords( c_vox', n, t, '' )' );
        end
    
        % print current row
        fprintf( 'Row %d out of %d finished!\n', i, nrow );
    
    end
    
    % save the result as .csv
    writetable( d, output_file )

end
