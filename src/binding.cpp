#include <nan.h>
#include <stdint.h>

class EncodeWorker : public Nan::AsyncWorker {
 public:
  EncodeWorker(
               Nan::Callback *end
  ) : Nan::AsyncWorker(end) {
  }

  ~EncodeWorker() {}

  void Execute () {
  }

 private:
  const uint8_t* context;
};

NAN_METHOD(create) {
}

NAN_METHOD(encode) {
}

NAN_MODULE_INIT(Init) {
  NAN_EXPORT(target, create); // Create an encoding context.
  NAN_EXPORT(target, encode); // Encode buffer or parity shards.
}

NODE_MODULE(binding, Init)

// S.D.G.
