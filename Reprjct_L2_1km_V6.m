% NOTE: NASA GSFC changed the structure and format of OceanColor SST data
% in late 2010.  Previously, 2030 control points were provided, one for
% each scan line (2030 scan control points, ie. cntl_pt_rows) but only 170
% pixel control points (ie. cntl_pt_cols) along the flight track.  This
% required the user to interpolate between the pixel control points.  A
% permanent change was implemented to have full lat/lon arrays on
% MODIS/Terra.  Evidently, the change was implemented in early 2010 for
% MODIS/Aqua.  (Ref: email from Bryan Franz dated 1/6/11 to Ryan).  The
% Reprjct code required updating to allow for a crude evaluation of pixel
% control point array size to determine which file structure type is used.


%%%%%%ENTER INPUT PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Infile directory & infile names
infiles=dir('*L2_LAC_SST.hdf');
proc_files=dir('*L2_SST_1km_jakob.mat');

%Outfile directory
outdirectory=['./'];

reg_tol=0.1; % tolerance for regional grid (%)
mic_tol=5;  % tolerance for micro grid or incremental grid  +/- value
radius_max=0.0468;       %maximum radius allowed
%   Parameters for Reprojected Region of Interest
LAT_S=68.5;
LAT_N=70.1;
LON_W=-53;
LON_E=-49.7;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Test for already processed images
list_all_files=char(zeros(size(infiles,1),14));
list_proc=char(zeros(size(proc_files,1),14));
list_unproc=char(zeros(size(list_all_files,1)-size(list_proc,1),14));
%read in list of processed files
for j=1:size(proc_files,1)
    list_proc(j,:)=proc_files(j).name(1:14);
end
%read in list of All files
for i=1:size(infiles,1)
    list_all_files(i,:)=infiles(i).name(1:14);
end
%Find number of UNprocessed files
kk=1;
for k=1:size(infiles,1)
    x=strmatch(list_all_files(k,:),list_proc(:,:));
    if (size(x,1)==0)
        list_unproc(kk,:)=list_all_files(k,:);
        kk=kk+1;
    end
end

% % Calculate number of steps & create grid for lat/lon reprojection - 1 km
% % spatial resolution
LAT_step=ceil((LAT_N-LAT_S)*111);
LON_step=ceil((LON_E-LON_W)*cosd(mean([LAT_S LAT_N]))*111);
new_lat=linspace(LAT_N,LAT_S,LAT_step);
new_lon=linspace(LON_W,LON_E,LON_step);
tic

for filenum=1:size(list_unproc,1)
    try
        clear numpixcntrlpts
        filename=[list_unproc(filenum,:),'.L2_LAC_SST.hdf'];
        fprintf(1,'Loading: %s\n',filename);
        %Read in data
        sst_info=hdfinfo(filename);
        numpixcntrlpts=sst_info.Vgroup(4).SDS(3).Dims.Size;
        sst_raw = single(hdfread(filename, 'sst', 'Index', {[1  1],[1  1],[2030  1354]}));
        qual_sst = hdfread(filename, 'qual_sst', 'Index', {[1  1],[1  1],[2030  1354]});
        l2_flags = hdfread(filename, 'l2_flags', 'Index', {[1  1],[1  1],[2030  1354]});
        longitude = hdfread(filename, 'longitude', 'Index', {[1  1],[1  1],[2030 numpixcntrlpts]});
        latitude = hdfread(filename, 'latitude', 'Index', {[1  1],[1  1],[2030  numpixcntrlpts]});
        cntl_pt_cols = hdfread(filename, 'cntl_pt_cols', 'Index', {[1],[1],[numpixcntrlpts]});
        cntl_pt_rows = hdfread(filename, 'cntl_pt_rows', 'Index', {[1],[1],[2030]});
        slope=sst_info.Vgroup(3).SDS(1).Attributes(2).Value;
        intercept=sst_info.Vgroup(3).SDS(1).Attributes(3).Value;
        SST=double(sst_raw)*slope+intercept;
        fprintf(1,'\nProcessing: %s,    file %d/%d\n',filename,filenum,size(list_unproc,1));
        
        if (numpixcntrlpts==170)
            %Interpolation routine to create lat/lon from L2 file
            lat_all=interp2(double(cntl_pt_cols),double(cntl_pt_rows),latitude,1:1354,[1:2030]');
            lon_all=interp2(double(cntl_pt_cols),double(cntl_pt_rows),longitude,1:1354,[1:2030]');
        else % lat/long data is in full array
            lat_all=latitude;
            lon_all=longitude;
        end
        
        %Calculate regional grid coordinates for region of interest
        dist_NW=sqrt((lat_all-(LAT_N*(1+reg_tol))).^2 + (cosd(lat_all).*(lon_all-(LON_W*(1+reg_tol))).^2));
        min_dist_NW=min(min(dist_NW));
        [i_NW,j_NW]=find(dist_NW==min_dist_NW);
        
        dist_NE=sqrt((lat_all-(LAT_N*(1+reg_tol))).^2 + (cosd(lat_all).*(lon_all-(LON_E*(1-reg_tol))).^2));
        min_dist_NE=min(min(dist_NE));
        [i_NE,j_NE]=find(dist_NE==min_dist_NE);
        
        dist_SW=sqrt((lat_all-(LAT_S*(1-reg_tol))).^2 + (cosd(lat_all).*(lon_all-(LON_W*(1+reg_tol))).^2));
        min_dist_SW=min(min(dist_SW));
        [i_SW,j_SW]=find(dist_SW==min_dist_SW);
        
        dist_SE=sqrt((lat_all-(LAT_S*(1-reg_tol))).^2 + (cosd(lat_all).*(lon_all-(LON_E*(1-reg_tol))).^2));
        min_dist_SE=min(min(dist_SE));
        [i_SE,j_SE]=find(dist_SE==min_dist_SE);
        
        i_4=[i_NW(1) i_NE(1) i_SE(1) i_SW(1)];j_4=[j_NW(1) j_NE(1) j_SE(1) j_SW(1)];
        min_ii=min(i_4);max_ii=max(i_4);min_jj=min(j_4);max_jj=max(j_4);
        
        %Initialize regional grids for radius calculations
        mini_lat=lat_all(min_ii:max_ii,min_jj:max_jj);
        mini_lon=lon_all(min_ii:max_ii,min_jj:max_jj);
        mini_sst=SST(min_ii:max_ii,min_jj:max_jj);
        mini_qual_sst=qual_sst(min_ii:max_ii,min_jj:max_jj);
        mini_l2_flags=l2_flags(min_ii:max_ii,min_jj:max_jj);
        
        %Initialize grids for newly projected data
        sst_proj=zeros(size(new_lat,2),size(new_lon,2));
        qual_sst_proj=int8(zeros(size(new_lat,2),size(new_lon,2)));
        l2_flags_proj=int32(zeros(size(new_lat,2),size(new_lon,2)));
        granule_lat=single(zeros(size(new_lat,2),size(new_lon,2)));
        granule_lon=single(zeros(size(new_lat,2),size(new_lon,2)));
        nn_radius=single(zeros(size(new_lat,2),size(new_lon,2)));
        
        %Initialize grids for micro projected data
        micro_sst=zeros(mic_tol*2+1,mic_tol*2+1);
        micro_qual_sst=int8(zeros(mic_tol*2+1,mic_tol*2+1));
        micro_l2_flags=int32(zeros(mic_tol*2+1,mic_tol*2+1));
        micro_granule_lat=single(zeros(mic_tol*2+1,mic_tol*2+1));
        micro_granule_lon=single(zeros(mic_tol*2+1,mic_tol*2+1));
        micro_nn_radius=single(zeros(mic_tol*2+1,mic_tol*2+1));
        
        ptlat=new_lat(1);
        ptlon=new_lon(1);
        dist=sqrt((mini_lat-ptlat).^2 + (cosd(mini_lat).*(mini_lon-ptlon)).^2);
        [min_vec,min_index_i]=min(dist);
        [min_val,min_index_j]=min(min_vec);
        min_j=min_index_j;
        min_i=min_index_i(min_index_j);
        loopcount=0;min_count=0;mic_count=0;scene_edgecount=0;
        lat_loop_time=zeros(1,size(new_lat,2));
        flag_scene_edge=0;
        
        for a=1:size(new_lat,2)
            move_on_to_new_line=0;
            for b=1:size(new_lon,2)
                loopcount=loopcount+1; %diagnostic counter
                ptlat=new_lat(a);
                ptlon=new_lon(b);
                
                mic_i=min_i; mic_j=min_j; mic_i_start= mic_i-mic_tol;
                mic_i_stop=mic_i+mic_tol; mic_j_start=mic_j-mic_tol; mic_j_stop=mic_j+mic_tol;
                
                if mic_i_start <= 0 || mic_i_start >= size(mini_lat,1) ...
                        || mic_i_stop >= size(mini_lat,1) ||mic_i_stop <= 0 ...
                        || mic_j_start <= 0 || mic_j_start >= size(mini_lat,2) ...
                        || mic_j_stop >= size(mini_lat,2) || mic_j_stop <= 0
                    
                    if mic_i_start <= 0
                        mic_i_start=1; mic_i_stop=2*mic_tol; edgeofgranule_flag=1;
                    end
                    if mic_i_start >= size(mini_lat,1)
                        mic_i_start= size(mini_lat,1)-mic_tol*2; mic_i_stop=size(mini_lat,1);
                        edgeofgranule_flag=1;
                    end
                    if mic_i_stop >= size(mini_lat,1)
                        mic_i_stop = size(mini_lat,1); mic_i_start = size(mini_lat,1) - mic_tol*2;
                        edgeofgranule_flag=1;
                    end
                    if mic_i_stop <= 0
                        mic_i_stop = 2*mic_tol; mic_i_start=1; edgeofgranule_flag=1;
                    end
                    
                    if mic_j_start <= 0
                        mic_j_start=1; mic_j_stop=mic_tol*2; edgeofgranule_flag=1;
                    end
                    if mic_j_start >= size(mini_lat,2)
                        mic_j_start=size(mini_lat,2); mic_j_stop=size(mini_lat,2)-2*mic_tol;
                        edgeofgranule_flag=1;
                    end
                    if mic_j_stop >= size(mini_lat,2)
                        mic_j_stop=size(mini_lat,2); mic_j_start=size(mini_lat,2)-2*mic_tol;
                        edgeofgranule_flag=1;
                    end
                    if mic_j_stop <= 0
                        mic_j_stop=2*mic_tol; mic_j_start=1; edgeofgranule_flag=1;
                    end
                end
                
                micro_sst=mini_sst(mic_i_start:mic_i_stop, mic_j_start:mic_j_stop);
                micro_lon=mini_lon(mic_i_start:mic_i_stop, mic_j_start:mic_j_stop);
                micro_lat=mini_lat(mic_i_start:mic_i_stop, mic_j_start:mic_j_stop);
                micro_qual_sst=mini_qual_sst(mic_i_start:mic_i_stop, mic_j_start:mic_j_stop);
                micro_l2_flags=mini_l2_flags(mic_i_start:mic_i_stop, mic_j_start:mic_j_stop);
                micro_granule_lat=mini_lat(mic_i_start:mic_i_stop, mic_j_start:mic_j_stop);
                micro_granule_lon=mini_lon(mic_i_start:mic_i_stop, mic_j_start:mic_j_stop);
                
                dist=sqrt((micro_lat-ptlat).^2 + (cosd(micro_lat).*(micro_lon-ptlon)).^2);
                [mic_vec,mic_index_i]=min(dist);
                [mic_val,mic_index_j]=min(mic_vec);
                mic_j=mic_index_j;
                mic_i=mic_index_i(mic_index_j);
                if (mic_val <= radius_max)
                    sst_proj(a,b)=micro_sst(mic_i,mic_j);
                    qual_sst_proj(a,b)=micro_qual_sst(mic_i,mic_j);
                    l2_flags_proj(a,b)=micro_l2_flags(mic_i,mic_j);
                    granule_lat(a,b)=micro_granule_lat(mic_i,mic_j);
                    granule_lon(a,b)=micro_granule_lon(mic_i,mic_j);
                    nn_radius(a,b)=mic_val;
                    min_i=min_i-mic_tol+mic_i;
                    min_j=min_j-mic_tol+mic_j;
                    mic_count=mic_count+1; %diagnostic counter
                else
                    dist=sqrt((mini_lat-ptlat).^2 + (cosd(mini_lat).*(mini_lon-ptlon)).^2);
                    [min_vec,min_index_i]=min(dist);
                    [min_val,min_index_j]=min(min_vec);
                    if (min_val <= radius_max)    % Loop to detect edge of scene
                        min_j=min_index_j;
                        min_i=min_index_i(min_index_j);
                        sst_proj(a,b)=mini_sst(min_i,min_j);
                        qual_sst_proj(a,b)=mini_qual_sst(min_i,min_j);
                        l2_flags_proj(a,b)=mini_l2_flags(min_i,min_j);
                        granule_lat(a,b)=mini_lat(min_i,min_j);
                        granule_lon(a,b)=mini_lon(min_i,min_j);
                        nn_radius(a,b)=min_val;
                        min_count=min_count+1; %diagnostic counter
                    else
                        sst_proj(a,b:size(new_lon,2))=999;
                        qual_sst_proj(a,b:size(new_lon,2))=999;
                        l2_flags_proj(a,b:size(new_lon,2))=999;
                        granule_lat(a,b:size(new_lon,2))=999;
                        granule_lon(a,b:size(new_lon,2))=999;
                        nn_radius(a,b:size(new_lon,2))=999;
                        scene_edgecount=scene_edgecount+1;
                        move_on_to_new_line=1;
                        flag_scene_edge=1;
                    end
                end
                if (move_on_to_new_line==1)
                    break;
                end
            end
        end
        
        clear ptlat ptlon dist min_vec min_index_i min_index_j min_val
        info1=['This data was created from L2 MODIS data and was reprojected using nearest neighbor criteria'];
        info2=['at 1-km spatial resolution.  The georef data contains the latitudes and longitudes for each of'];
        info3=['the newly projected SST data.  The granule_lat/granule_lon contains interpolated lat/long from'];
        info4=['the original granule. The flag_scene_edge identifies reprocessed scenes with granule edge or '];
        info5=['missing data present.  A zero represents a normal scene or not missing data. A one represents a '];
        info6=['granule edge of missing data is present'];
        info=[info1,sprintf('\n'),info2,sprintf('\n'),info3,sprintf('\n'),info4,sprintf('\n'),info5,sprintf('\n'),info6];
        version=['This data processed using the Reprjct_L2_1km_V5a_jakob.m script'];
        
        L2_data.filename=filename;
        L2_data.sst=sst_proj;
        L2_data.georef.lat=new_lat;
        L2_data.georef.lon=new_lon;
        L2_data.nn_radius=nn_radius;
        L2_data.qual_sst=qual_sst_proj;
        L2_data.l2_flags=l2_flags_proj;
        L2_data.granule_lat=granule_lat;
        L2_data.granule_lon=granule_lon;
        L2_data.info=info;
        L2_data.version=version;
        L2_data.flag_scene_edge=flag_scene_edge;
        L2_data.date_created=date;
        
        outfilename=[outdirectory,filename(1:18),sprintf('%c'),'SST_1km_jakob.mat'];
        fprintf(1,'Saving %s\n',outfilename);
        save(outfilename,'L2_data')
        toc
        clear L2_data filename outfilename dist* max* lat* lon* info* gran* sst* qual* l2* i* SST cntl*
        
    catch
        fprintf(2,'\n%s did not get processed.  Granule mismatch.\n\n', filename)
        error_name=filename(1:25);
        OUTOFBOUNDS=['_OUTOFBOUNDS'];
        suffix=['_hdf.txt'];
        errorfilename=[error_name,OUTOFBOUNDS,suffix];
        error_msg=['!!Granule did not process properly.  Likely a granule mismatch!!'];
        dlmwrite(errorfilename, error_msg, 'delimiter', '')
        continue
    end
end

fprintf(2,'\n\n...........All images in directory have been processed\n')
