function [intensity] = point_sources_intensity(pressure,velocity)
%POINT_SOURCES_INTENSITY Calculate the intensity field of point sources.

intensity = (1/2)*pressure.*conj(velocity);
end


