```bash
# AddressSanitizer
docker build -f Dockerfile.biobambam2-debug --tag biobambam2-debug-asan --build-arg CCFLAGS=-fsanitize=address --build-arg CXXFLAGS=-fsanitize=address .

# ThreadSanitizer
docker build -f Dockerfile.biobambam2-debug --tag biobambam2-debug-tsan --build-arg CCFLAGS=-fsanitize=thread --build-arg CXXFLAGS=-fsanitize=thread .
```
