function [start_light, end_light,dataStruct, timeAxisStruct]= resize3 (Light_trace, dataStruct, samplingInts, timeAxisStruct, percent)
fields= fieldnames(dataStruct);
s=length(fields);
si=dataStruct.(fields{1});
siz=size(si);
num_ch=siz(2);

start_light=[];
end_light=[];

%intracellular 1 and 2;

%Light_trace=4;


for i=1:s
    si=dataStruct.(fields{i});
    sa(i)=size(si,1);
end
mi_s=min(sa);
for i=1:s
    si=dataStruct.(fields{i});
    si=si(1:mi_s,:);
    dataStruct.(fields{i})=si;
    si_flip=fliplr(si');
    si_flip=si_flip';
    duration(i)=length(si(:,1));
    d_light=diff(si(:,Light_trace));
    d_light_flip=diff(si_flip(:, Light_trace));
    m=find(d_light>1);
    m=min(m);%max(d_light);
    if ~isempty(m)
    k1(i)=m;
    else
        k1=[];
    end
    
    m_flip=find(d_light_flip>1);

    m=min(m_flip);%min(d_light);
   if ~isempty(m)
    %k=find(m==d_light);
    k2(i)=duration(i)-m;
   else
       k2=[];
   end
    
   
end

if ~isempty(k1)
start_light=min(k1);
end
if ~isempty(k1) && ~isempty(k2)
end_light=max(k2);
end

if ~isempty(k1)
light_duration=(((k2-k1)*samplingInts(1)));
start_light=(percent*start_light/100);
temp1=floor(start_light);
start_light=(start_light*samplingInts(1));
temp2=floor((duration(1)-end_light)*percent/100);
%temp3=(round(end_light*samplingInts(1)/10))*10/samplingInts(1);
end_light=start_light+light_duration;
%start_light=end_light-light_duration;

%end_light_f=end_light/samplingInts(1);
%start_light_f=start_light/samplingInts(1);


for i=1:s
    %back_keep=duration(i)-end_light;
  %  back_keep=ceil(50/samplingInts(1));%50
si=dataStruct.(fields{i});

    %si(k1(i):k2(i),Light_trace)=1;
    %si(temp1:k1, Light_trace)=0;
   % si(k2(i):temp2, Light_trace)=0;
si_cut=si(k1(i)-temp1:k2(i)+temp2,:);%50

l=length(si_cut);
%si_cut=resample(si_cut, 800000, l);
dataStruct.(fields{i})=si_cut;
%ti=timeAxisStruct.(fields{i});

%ti_cut=ti(k1(i)-temp1:k2(i)+temp2);%0:samplingInts(i):(l*samplingInts(i));
ti=linspace(0,l*samplingInts(i), size(si_cut,1) );
%ti=resample(ti, 800000, l);
timeAxisStruct.(fields{i})=ti;
end
end
end
