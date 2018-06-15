function taskStartup(task)

%% Function to setup the rng for the different workers based on a substream system
% Based on discussion in https://www.mathworks.com/matlabcentral/newsreader/view_thread/327817

% Get task ID, (parent) job ID and job creation time
taskID = task.ID;
job = task.Parent;
jobID = job.ID;
job_creation_time=job.CreateTime;

% Convert the job creation time to seconds
job_creation_time=job_creation_time([1:20,25:28]); %remove the time zone
job_creation_time=round(datenum(job_creation_time,'ddd mmm dd HH:MM:SS yyyy')*86400); %convert (units: seconds)

% Create shift indices using the job creation time:
shift_index_stream=job_creation_time; %shift index for stream
shift_index_substream=round(job_creation_time/1000); %shift index for substream

%Now:
%1) Create a large number of independent streams;
%2) Select the stream using shift_index_stream and jobID
%3) For this particular task use a substream identified by
% shift_index_substream and taskID
%

NS=2^63; %number of multiple independent streams
s = RandStream.create('mrg32k3a', 'NumStreams', NS, 'Stream', shift_index_stream+jobID); %note: there is no shuffle
s.Substream = shift_index_substream+taskID;
RandStream.setGlobalStream(s);

% This process does not use a new seed, which could cause issues with
% independence.
% Instead, it relies on a unique identifier based on the task and job id,
% plus an rng element based on the job creation time, in order to select
% out a particular substream that was generated from the right
% distribution.
end