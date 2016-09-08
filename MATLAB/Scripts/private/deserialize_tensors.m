function [C, T] = deserialize_tensors(record)
  C = getArrayFromByteStream(record.c);
  T = getArrayFromByteStream(record.t);
end
