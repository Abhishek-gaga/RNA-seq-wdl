IMAGE_NAME = bacterial-rnaseq-umi:1.0.0
VERSION = 1.0.0
TAG1 = Abhishek-gaga/$(IMAGE_NAME):$(VERSION)
TAG2 = Abhishek-gaga/$(IMAGE_NAME):latest

all: | build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
