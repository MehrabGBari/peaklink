function runNotification(msg)


if exist('notification.p', 'file') ~= 0;
    notification(msg)
else
    disp('notification not performed.');
end